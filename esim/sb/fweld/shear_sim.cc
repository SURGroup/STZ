#include <cstring>

#include "common.hh"
#include "shear_sim.hh"
#include "mat.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *		      vertical directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] mu_ the elastic shear modulus.
 * \param[in] K_ the elastic bulk modulus.
 * \param[in] viscosity_ the viscosity.
 * \param[in] t_scale_ a time scale associated with the plastic deformation.
 * \param[in] ub_ the boundary velocity.
 * \param[in] tmult_ a multiplier to apply to the default timestep size.
 * \param[in] adapt_fac_ a factor used in the adaptive timestepping of
 *			 plasticity.
 * \param[in] filename_ the filename of the output directory. */
shear_sim::shear_sim(const int m_,const int n_,const double ax_,const double bx_,
		const double ay_,const double by_,const double mu_,
		const double K_,const double viscosity_,const double t_scale_,
		const double tmult_,const double adapt_fac_,const double u_bdry_,
		stz_dynamics *stz_,const char *filename_)
	: m(m_), n(n_), mn(m_*n_), ax(ax_), ay(ay_), bx(bx_), by(by_),
	dx((bx_-ax_)/m_), dy((by_-ay_)/n_), xsp(1/dx), ysp(1/dy), mu(mu_),
	mu_inv(1/mu_), K(K_), viscosity(viscosity_), t_scale(t_scale_), tmult(tmult_),
	adapt_fac(adapt_fac_), u_bdry(u_bdry_), bcs(4), stz(stz_), filename(filename_),
	fm((new c_field[mn+2*m])+m), vel(new vec[mn+m]), src(new vec[mn+m]), qsm(*this),
	mg(qsm,src,vel), buf(new float[m>=63?m+2:64]) {

	// Set STZ related parameters
	TZ=stz->TZ;
	chi_inf=stz->chi_inf;
}

/** The class destructor frees the dynamically allocated memory. */
shear_sim::~shear_sim() {
	delete [] buf;
	delete [] src;
	delete [] vel;
	delete [] (fm-m);
}

/** Initializes the simulation fields.
 * \param[in] chi_case the type of setup to use for chi.
 * \param[in] tem_base the base value of effective temperature (in Kelvin) to
 *		       use.
 * \param[in] tem_delta the delta value of effective temperature (in Kelvin) to
 *			use. */
void shear_sim::init_fields(int chi_case,double tem_base,double tem_delta) {
	const double pi=3.1415926535897932384626433832795;
	double chi0=tem_base/TZ,chi1=tem_delta/TZ,xdis=bx-ax;

	// Initialize table for sine wave case
	double sx[16],sy[16];
	if(chi_case==2) {
		for(int i=0;i<16;i++) {
			sx[i]=ax+xdis*(i+0.5)*0.0625;
			sy[i]=0.2*sin(pi*0.0625*(2*i+1));
		}
	}
#pragma omp parallel for
	for(int j=-1;j<n+1;j++) {
		double Y=ay+j*dy,Ys=Y+0.5*dy;
		c_field *fp=fm+m*j;
		for(int i=0;i<m;i++,fp++) {
			c_field &f=*fp;
			double X=ax+i*dx,Xs=X+0.5*dx;

			// Set stresses and velocities to zero
			f.u=f.v=0;
			f.p=f.q=f.s=f.tau=0;

			// Zeroth case: a small Gaussian blip in chi at the
			// origin
			if(chi_case==0) f.chi=chi0+chi1*exp(-200*(Xs*Xs+Ys*Ys));
			else if(chi_case==1) {

				// First case: a rotated line of higher chi
				const double csa=cos(30*pi/180.),sia=sin(30*pi/180.);
				double Xr=csa*Xs+sia*Ys,Yr=-sia*Xs+csa*Ys;
				Xr=Xr<-1?Xr+1:(Xr>1?Xr-1:0);
				f.chi=chi0+chi1*exp(-200*(Xr*Xr+Yr*Yr));
			} else {

				// Second case: a sequence of Gaussian blips in
				// a sine wave pattern
				double R=1e100,Rt,Xr,Yr;
				for(int k=0;k<16;k++) {
					Xr=sx[k]-Xs;
					if(Xr<-xdis*0.5) Xr+=xdis;
					if(Xr>xdis*0.5) Xr-=xdis;
					Yr=sy[k]-Ys;
					Rt=Xr*Xr+Yr*Yr;
					if(Rt<R) R=Rt;
				}
				f.chi=chi0+chi1*exp(-200*R);
			}

			// Set reference map field to match initial coordinate
			f.X=X;f.Y=Y;
		}
	}
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output.
 * \param[in] t_start the simulation time to start from.
 * \param[in] t_end the simulation time to end at.
 * \param[in] frames the number of frames to save. */
void shear_sim::solve(double t_start,double t_end,int frames) {
	time=t_start;
	double time_interval=(t_end-t_start)/frames,target_time,t0,t1,t2;
	int l;
	const double dt=dx*dx*tmult;

	// Output the initial fields and record initial time
	write_files(0);
	puts("# Output frame 0");
	t0=wtime();

	for(int k=1;k<=frames;k++) {

		// Compute the target time to the next output frame
		target_time=t_start+time_interval*k;l=1;

		// Carry out simulation step using the regular timestep until
		// within range of the target time
		while(time+dt*(1+1e-8)<target_time) {l++;step_forward(dt);}

		// Carry out a final simulation step, using exactly the right
		// time step to reach the target time
		step_forward(target_time-time);

		// Output the fields
		t1=wtime();
		write_files(k);

		// Print diagnostic information
		t2=wtime();
		printf("# Output frame %d [%d, %.8g s, %.8g s]\n",k,l,t1-t0,t2-t1);
		t0=t2;
	}
}

/** Carries out the simulation for a specified time interval using the
 * quasi-static simulation method, periodically saving the output.
 * \param[in] t_start the simulation time to start from.
 * \param[in] t_end the simulation time to end at.
 * \param[in] frames the number of frames to save.
 * \param[in] steps the number of quasi-static steps to take per frame. */
void shear_sim::solve_quasistatic(double t_start,double t_end,int frames,int steps) {
	time=t_start;
	double time_interval=(t_end-t_start)/frames,t0,t1,t2;
	double isteps=1.0/steps;

	// Output the initial fields
	write_files(0);
	puts("# Output frame 0");
	t0=wtime();

	for(int k=1;k<=frames;k++) {

		for(int j=0;j<steps;j++) step_forward_quasistatic(time_interval*isteps);
		t1=wtime();

		// Output the fields
		write_files(k);

		// Print diagnositic information
		t2=wtime();
		printf("# Output frame %d [%d, %.8g s, %.8g s]\n",k,steps,t1-t0,t2-t1);
		t0=t2;
	}
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void shear_sim::step_forward(double dt) {
	const double third=1./3;
	const double lx=bx-ax;
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt;
	double hxxv=2*xsp*hx*viscosity,hyyv=2*ysp*hy*viscosity,hcv=2*(hxxv+hyyv);

	// Update the simulation time
	time+=dt;
//	printf("%g\n",time);

#pragma omp parallel for
	for(int j=0;j<n;j++) for(int i=0;i<m;i++) {
		st_field f1,f2;

		// Create references to the fields in the neighboring gridpoints,
		// taking into account the periodicity
		c_field *fp=fm+(m*j+i),*flp=fp+(i==0?m-1:-1),*frp=fp+(i==m-1?1-m:1);
		c_field &f=*fp,&fl=*flp,&fr=*frp,&fd=fp[-m],&fu=fp[m];
		double uc,vc,ux,vx,uy,vy,uvisc,vvisc,two_omega,ddev;

		// Unless this is the zeroth row, calculate the update to the
		// reference map and velocity
		if(j>0) {
			double Xx,Xy,Yx,Yy,fx,fy;
			uc=f.u;vc=f.v;

			// Calculate ENO derivatives of regular fields
			uc>0?rmv_eno2(ux,vx,Xx,Yx,hx,fr,f,fl,fp[i<=1?m-2:-2],
				      i==m-1?lx:0,i==0?-lx:0,i<=1?-lx:0)
			    :rmv_eno2(ux,vx,Xx,Yx,-hx,fl,f,fr,fp[i>=m-2?2-m:2],
				      i==0?-lx:0,i==m-1?lx:0,i>=m-2?lx:0);
			vc>0?(j>1?rmv_eno2(uy,vy,Xy,Yy,hy,fu,f,fd,fp[-2*m])
				 :rmv_one_sided(uy,vy,Xy,Yy,hy,f,fd))
			    :(j<n-1?rmv_eno2(uy,vy,Xy,Yy,-hy,fd,f,fu,fp[2*m])
				   :rmv_one_sided(uy,vy,Xy,Yy,-hy,f,fu));

			// Calculate net force due to stress imbalance
			net_force(fx,fy,-hx,-hy,flp,fp);

			// Calculate viscosity
			uvisc=hxxv*(fr.u+fl.u)+hyyv*(fu.u+fd.u)-hcv*uc;
			vvisc=hxxv*(fr.v+fl.v)+hyyv*(fu.v+fd.v)-hcv*vc;

			// Apply the update to velocity
			f.cu=-uc*ux-vc*uy+mu_inv*(fx+uvisc);
			f.cv=-uc*vx-vc*vy+mu_inv*(fy+vvisc);

			// Calculate update to reference map fields
			f.cX=-uc*Xx-vc*Xy;
			f.cY=-uc*Yx-vc*Yy;
		}

		// Calculate deformation rate, staggered velocity, and spin
		uc=0.25*(f.u+fr.u+fu.u+frp[m].u);
		vc=0.25*(f.v+fr.v+fu.v+frp[m].v);
		ux=hx*(-f.u+fr.u-fu.u+frp[m].u);
		uy=hy*(-f.u-fr.u+fu.u+frp[m].u);
		vx=hx*(-f.v+fr.v-fu.v+frp[m].v);
		vy=hy*(-f.v-fr.v+fu.v+frp[m].v);
		two_omega=vx-uy;

		// Calculate ENO derivatives of staggered fields
		f1=uc>0?st_eno2(hx,fr,f,fl,fp[i<=1?m-2:-2]):st_eno2(-hx,fl,f,fr,fp[i>=m-2?2-m:2]);
		f2=vc>0?(j>0?st_eno2(hy,fu,f,fd,fp[-2*m]):st_one_sided(hy,f,fd))
		       :(j<n-1?st_eno2(-hy,fd,f,fu,fp[2*m]):st_one_sided(-hy,f,fu));

		// Calculate the adaptive effective temperature update
		f.cchi=f.chi;
		ddev=adaptive_plastic_term(f.dev(),f.cchi,dt);
		f.cchi-=uc*f1.chi+vc*f2.chi;

		// Calculate updates to stress
		f.cp=-uc*f1.p-vc*f2.p-K*(ux+vy);
		f.cq=-uc*f1.q-vc*f2.q-mu*third*(ux+vy)-ddev*f.q;
		f.cs=-uc*f1.s-vc*f2.s-two_omega*f.tau+mu*(ux-vy)-ddev*f.s;
		f.ctau=-uc*f1.tau-vc*f2.tau+two_omega*f.s+mu*(vx+uy)-ddev*f.tau;
	}

	// Apply the updates to the fields
#pragma omp parallel for
	for(int j=0;j<n;j++) {
		if(j>0) {for(int i=0;i<m;i++) fm[i+m*j].update();}
		else for(int i=0;i<m;i++) fm[i+m*j].update_staggered();
	}

	// Set the top and bottom boundaries
	set_boundaries();
}

/** Steps the simulation fields forward using the quasistatic update procedure.
 * \param[in] dt the time step to use. */
void shear_sim::step_forward_quasistatic(double dt) {

	// Update the simulation time
	time+=dt;

	// Carry out the intermediate step that considers the advective and
	// plasiticty terms
	advection_step(dt);

	// Apply the projection step to ensure quasistaticity
	projection_step(dt);
}

/** Carries out an intermediate step forward, by evaluating the advective and
 * plasticity terms in the equations only.
 * \param[in] dt the time step to use. */
void shear_sim::advection_step(double dt) {
	const double lx=bx-ax;
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt;

#pragma omp parallel for
	for(int j=0;j<n;j++) for(int i=0;i<m;i++) {
		st_field f1,f2;

		// Create references to the fields in the neighboring gridpoints,
		// taking into account the periodicity
		c_field *fp=fm+(i+m*j),*flp=fp+(i==0?m-1:-1),*frp=fp+(i==m-1?1-m:1);
		c_field &f=*fp,&fl=*flp,&fr=*frp,&fd=fp[-m],&fu=fp[m];
		double uc,vc,two_omega,ddev;

		// Calculate update to reference map fields
		if(j>0) {
			double Xx,Xy,Yx,Yy;
			uc=f.u;vc=f.v;

			// Calculate ENO derivatives of the reference map
			uc>0?rm_eno2(Xx,Yx,hx,fr,f,fl,fp[i<=1?m-2:-2],
				     i==m-1?lx:0,i==0?-lx:0,i<=1?-lx:0)
			    :rm_eno2(Xx,Yx,-hx,fl,f,fr,fp[i>=m-2?2-m:2],
				     i==0?-lx:0,i==m-1?lx:0,i>=m-2?lx:0);
			vc>0?(j>1?rm_eno2(Xy,Yy,hy,fu,f,fd,fp[-2*m])
				 :rm_one_sided(Xy,Yy,hy,f,fd))
			    :(j<n-1?rm_eno2(Xy,Yy,-hy,fd,f,fu,fp[2*m])
				   :rm_one_sided(Xy,Yy,-hy,f,fu));

			// Calculate update to reference map fields
			f.cX=-uc*Xx-vc*Xy;
			f.cY=-uc*Yx-vc*Yy;
		}

		// Calculate deformation rate, staggered velocity, and spin
		uc=0.25*(f.u+fr.u+fu.u+frp[m].u);
		vc=0.25*(f.v+fr.v+fu.v+frp[m].v);
		two_omega=hx*(-f.v+fr.v-fp[m].v+frp[m].v)
			 -hy*(-f.u-fr.u+fp[m].u+frp[m].u);

		// Calculate ENO derivatives of staggered fields
		f1=uc>0?st_eno2(hx,fr,f,fl,fp[i<=1?m-2:-2]):st_eno2(-hx,fl,f,fr,fp[i>=m-2?2-m:2]);
		f2=vc>0?(j>0?st_eno2(hy,fu,f,fd,fp[-2*m]):st_one_sided(hy,f,fd))
		       :(j<n-1?st_eno2(-hy,fd,f,fu,fp[2*m]):st_one_sided(-hy,f,fu));

		// Calculate the adaptive effective temperature update
		f.cchi=f.chi;
		ddev=adaptive_plastic_term(f.dev(),f.cchi,dt);
		f.cchi-=uc*f1.chi+vc*f2.chi;

		// Calculate updates to stresses
		f.cp=-uc*f1.p-vc*f2.p;
		f.cq=-uc*f1.q-vc*f2.q-ddev*f.q;
		f.cs=-uc*f1.s-vc*f2.s-two_omega*f.tau-ddev*f.s;
		f.ctau=-uc*f1.tau-vc*f2.tau+two_omega*f.s-ddev*f.tau;
	}

	// Apply the updates to the fields
#pragma omp parallel for
	for(int j=0;j<n;j++) {
		if(j>0) {for(int i=0;i<m;i++) fm[i+m*j].update_others();}
		else for(int i=0;i<m;i++) fm[i+m*j].update_staggered();
	}
}

/** Carries out a projection step of the stresses using a multigrid solve,
 * to enforce the quasistaticity constraint.
 * \param[in] dt the timestep to use. */
void shear_sim::projection_step(double dt) {
	const double third=1/3.0;
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt,xfac=0.5/dt*xsp,yfac=0.5/dt*ysp;
	int i,j;

	// Initialize constants in the algebraic system
	qsm.init(viscosity/dt);

	// Set lower boundary condition
	double xbc=bcs*qsm.cxx;
	for(i=0;i<m;i++) {src[i].x=-xbc*bdry_vel();src[i].y=0;}

	// Calculate the source term for the projection step
#pragma omp parallel for
	for(j=1;j<n;j++) for(int ij=m*j;ij<(j+1)*m;ij++)
		net_force(src[ij].x,src[ij].y,xfac,yfac,fm+ij+(ij==m*j?m-1:-1),fm+ij);

	// Set upper boundary condition
	for(i=mn;i<mn+m;i++) {src[i].x=xbc*bdry_vel();src[i].y=0;}

	// Set up the algebraic systems in the grid hierarchy and carry out the
	// multigrid solve
	mg.setup();
	mg.solve_v_cycle();

	// Apply the correction to the stress tensor, based upon the computed
	// velocity. In addition, set up a representation of the velocity on
	// the primary grid
#pragma omp parallel for
	for(j=0;j<n;j++) {
		for(int i=0;i<m;i++) {
			c_field &f=fm[i+m*j];
			vec *vp=vel+(i+m*j),*vo=vp+(i==m-1?1-m:1);
			f.u=vp->x;
			f.v=vp->y;
			double ux=hx*(-vp->x+vo->x-vp[m].x+vo[m].x),
			       uy=hy*(-vp->x-vo->x+vp[m].x+vo[m].x),
			       vx=hx*(-vp->y+vo->y-vp[m].y+vo[m].y),
			       vy=hy*(-vp->y-vo->y+vp[m].y+vo[m].y);
			f.p+=-K*(ux+vy);
			f.q+=-mu*third*(ux+vy);
			f.s+=mu*(ux-vy);
			f.tau+=mu*(uy+vx);
		}
	}

	// Set the boundaries again
	set_boundaries();
}

/** Calculates the net force at a regular grid point using the stress tensors
 * on the four neighboring points of the staggered grids.
 * \param[out] (fx,fy) the components of force.
 * \param[in] (xf,yf) prefactors to apply to the x and y derivatives that go into the computation of force.
 * \param[in] fp0 a pointer to the grid point to the left of the one to
 * 		  consider.
 * \param[in] fp1 a pointer to the grid point to consider. */
inline void shear_sim::net_force(double &fx,double &fy,double xf,double yf,c_field *fp2,c_field *fp3) {
	c_field &f0=fp2[-m],&f1=fp3[-m],&f2=*fp2,&f3=*fp3;
	fx=(-f0.p+f1.p-f2.p+f3.p-f0.q+f1.q-f2.q+f3.q+f0.s-f1.s+f2.s-f3.s)*xf
	  +(f0.tau+f1.tau-f2.tau-f3.tau)*yf;
	fy=(f0.tau-f1.tau+f2.tau-f3.tau)*xf
	  +(-f0.p-f1.p+f2.p+f3.p-f0.q-f1.q+f2.q+f3.q-f0.s-f1.s+f2.s+f3.s)*yf;
}

/** Evaluates the squared L2 norm between the simulation fields and the
 * simulation fields in another instance of the shear_sim class.
 * \param[in] ss another shear_sim class.
 * \param[in] l2 a pointer to a four-entry array in which to store the output.
 */
void shear_sim::l2_comparison(shear_sim &ss,double *l2) {
	*l2=l2[1]=l2[2]=l2[3]=0;

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		c_field *fp=fm+m*j,*fo=ss.fm+m*j,*fe=fp+m;
		double c0=0,c1=0,c2=0,c3=0;
		while(fp<fe) {
			double du=fp->u-fo->u,dv=fp->v-fo->v;
			double dp=fp->p-fo->p,dq=fp->q-fo->q;
			double ds=fp->s-fo->s,dtau=fp->tau-fo->tau;
			double dchi=fp->chi-fo->chi;
			double dX=fp->X-fo->X,dY=fp->Y-fo->Y;
			c0+=du*du+dv*dv;
			c1+=2*(dtau*dtau+ds*ds)+3*dp*dp+6*dq*dq;
			c2+=dchi*dchi;
			c3+=dX*dX+dY*dY;
			fp++;fo++;
		}

		// Scale the zeroth line of the regular fields according to the
		// trapezoidal rule
		if(j==0) {c0*=0.5;c3*=0.5;}

		// Add the totals for this line to the accumulators
#pragma omp critical
		{*l2+=c0;l2[1]+=c1;l2[2]+=c2;l2[3]+=c3;}
	}

	// Deal with the extra line for the regular fields
	c_field *fp=fm+mn,*fo=ss.fm+mn,*fe=fp+m;
	while(fp<fe) {
		double du=fp->u-fo->u,dv=fp->v-fo->v;
		double dX=fp->X-fo->X,dY=fp->Y-fo->Y;
		*l2+=0.5*(du*du+dv*dv);
		l2[3]+=0.5*(dX*dX+dY*dY);
		fp++;fo++;
	}

	// Normalize by the grid size
	double fac=dx*dy;
	*l2*=fac;l2[1]*=fac;l2[2]*=fac;l2[3]*=fac;
}

/** Adaptively calculates the effect of plasticity on the shear stress
 * components and the effective temperature, by dividing up the timestep
 * into smaller intervals to ensure that the change in shear stress remains
 * small.
 * \param[in] sbar the magnitude of the shear stress at a grid point.
 * \param[in,out] chiv the effective temperature at the grid point.
 * \param[in] dt the timetep.
 * \return The factor to scale the shear stresses by. */
double shear_sim::adaptive_plastic_term(double sbar,double &chiv,double dt) {
	const double small_dev_cutoff=1e-12;
	const double theta=0;
	const bool no_chidot=false;
	bool adapt=true;
	double dtl=dt,osbar=sbar,dplas,dchi1,dchi2,adt;

	// To avoid division-by-zero errors, check for the case when the
	// deviatoric stress is small, and just return with no changes
	if(sbar<small_dev_cutoff) return 0;
	do {

		// Calculate the plastic deformation
		dplas=2*t_scale*mu*stz->Dplastic(sbar,chiv,dchi1,dchi2);

		// If the plastic deformation is large, then adaptively
		// choose the timestep so as not to make a large alteration in
		// the deviatoric stress magnitude
		if(abs(dtl*dplas)>adapt_fac) {
			adt=adapt_fac/abs(dplas);
			dtl-=adt;
		} else {
			adapt=false;
			adt=dtl;
		}
		sbar-=adt*dplas;

		// Calculate the change in chi if needed
		if(!no_chidot) chiv+=(dchi1*(chi_inf-chiv)+dchi2*(theta-chiv))*t_scale*adt;
	} while(adapt);

	// Calculate the quantity dt*2*mu*Dpl/sbar, which will be used in the
	// finite-difference update
	return 1-sbar/osbar;
}

/** Sets the fields in the top and bottom boundaries according to the boundary
 * conditions. */
void shear_sim::set_boundaries() {

	// Set bottom boundary according to the boundary conditions
	for(int i=0;i<m;i++) {
		c_field *fp=fm+i,&f=fp[-m],&fu=*fp,&fuu=fp[m];
		fu.u=-bdry_vel();
		fu.v=0;
		fu.X=ax+i*dx+bdry_pos();
		fu.Y=ay;
		f.p=2*fu.p-fuu.p;
		f.q=2*fu.q-fuu.q;
		f.s=2*fu.s-fuu.s;
		f.tau=2*fu.tau-fuu.tau;
		f.chi=2*fu.chi-fuu.chi;
	}

	// Set top boundary according to the boundary conditions
	for(int i=0;i<m;i++) {
		c_field *fp=fm+(mn+i),&f=*fp,&fd=fp[-m],&fdd=fp[-2*m];
		f.u=bdry_vel();
		f.v=0;
		f.X=ax+i*dx-bdry_pos();
		f.Y=by;
		f.p=2*fd.p-fdd.p;
		f.q=2*fd.q-fdd.q;
		f.s=2*fd.s-fdd.s;
		f.tau=2*fd.tau-fdd.tau;
		f.chi=2*fd.chi-fdd.chi;
	}
}

/** Calculates one-sided derivatives of the fields using first-order
 * upwinding.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f1,f2) the fields to compute the derivative with.
 * \return The computed derivative. */
inline void shear_sim::rmv_one_sided(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f1,c_field &f2) {
	ud=2*hs*(f1.u-f2.u);
	vd=2*hs*(f1.v-f2.v);
	rm_one_sided(Xd,Yd,hs,f1,f2);
}

/** Calculates one-sided derivatives of the fields using first-order
 * upwinding.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f1,f2) the fields to compute the derivative with.
 * \return The computed derivative. */
inline void shear_sim::rm_one_sided(double &Xd,double &Yd,double hs,c_field &f1,c_field &f2) {
	hs*=2;
	Xd=hs*(f1.X-f2.X);
	Yd=hs*(f1.Y-f2.Y);
}

/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme, applying the shift to the X terms.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \param[in] (X0,X2,X3) corresponding shifts to apply to the X terms.
 * \return The computed derivative. */
inline void shear_sim::rmv_eno2(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3,double X0,double X2,double X3) {
	ud=hs*eno2(f0.u,f1.u,f2.u,f3.u);
	vd=hs*eno2(f0.v,f1.v,f2.v,f3.v);
	rm_eno2(Xd,Yd,hs,f0,f1,f2,f3,X0,X2,X3);
}
/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme, applying the shift to the X terms.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \param[in] (X0,X2,X3) corresponding shifts to apply to the X terms.
 * \return The computed derivative. */
inline void shear_sim::rm_eno2(double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3,double X0,double X2,double X3) {
	Xd=hs*eno2(f0.X+X0,f1.X,f2.X+X2,f3.X+X3);
	Yd=hs*eno2(f0.Y,f1.Y,f2.Y,f3.Y);
}

/** Calculates one-sided derivatives of the fields using first-order
 * upwinding.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f1,f2) the fields to compute the derivative with.
 * \return The computed derivative. */
inline st_field shear_sim::st_one_sided(double hs,c_field &f1,c_field &f2) {
	hs*=2;
	return st_field(hs*(f1.p-f2.p),hs*(f1.q-f2.q),hs*(f1.s-f2.s),hs*(f1.tau-f2.tau),
		hs*(f1.chi-f2.chi));
}

/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme, applying the shift to the X terms.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \return The computed derivative. */
inline st_field shear_sim::st_eno2(double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3) {
	return st_field(hs*eno2(f0.p,f1.p,f2.p,f3.p),hs*eno2(f0.q,f1.q,f2.q,f3.q),
		hs*eno2(f0.s,f1.s,f2.s,f3.s),hs*eno2(f0.tau,f1.tau,f2.tau,f3.tau),
		hs*eno2(f0.chi,f1.chi,f2.chi,f3.chi));
}

/** Calculates the ENO derivative using a sequence of values at
 * four gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double shear_sim::eno2(double p0,double p1,double p2,double p3) {
	return abs(p0-2*p1+p2)>abs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void shear_sim::write_files(int k) {
	const int fflags=4|128|256;//512|1024;
	if(fflags&1) output("u",0,k);
	if(fflags&2) output("v",1,k);
	if(fflags&4) output("p",2,k);
	if(fflags&8) output("q",3,k);
	if(fflags&16) output("s",4,k);
	if(fflags&32) output("tau",5,k);
	if(fflags&64) output("chi",6,k);
	if(fflags&128) output("tem",7,k);
	if(fflags&256) output("dev",8,k);
	if(fflags&512) output("X",9,k);
	if(fflags&1024) output("Y",10,k);
	if(fflags&2048) {
		double Q,max_qs;
		qs_measure(Q,max_qs);
		printf("QS info %g %.12g %.12g\n",time,Q,max_qs);
		output("proj",11,k);
	}
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void shear_sim::output(const char *prefix,const int mode,const int sn) {

	// Determine whether to output a staggered field or not
	bool st=mode>=2&&mode<=8;
	int l=st?m:m+1;

	// Assemble the output filename and open the output file
	char *bufc=((char*) buf);
	sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
	FILE *outf=safe_fopen(bufc,"wb");

	// Output the first line of the file
	int i,j;
	float *bp=buf+1,*be=buf+m+1;
	*buf=l;
	for(i=0;i<l;i++) *(bp++)=ax+(st?(i+0.5)*dx:i*dx);
	fwrite(buf,sizeof(float),l+1,outf);

	// Output the field values to the file
	c_field *fp=fm;
	for(j=0;j<(st?n:n+1);j++) {
		*buf=ay+(st?(j+0.5)*dy:j*dy);bp=buf+1;
		switch(mode) {
			case 0: while(bp<be) *(bp++)=(fp++)->u;break;
			case 1: while(bp<be) *(bp++)=(fp++)->v;break;
			case 2: while(bp<be) *(bp++)=(fp++)->p;break;
			case 3: while(bp<be) *(bp++)=(fp++)->q;break;
			case 4: while(bp<be) *(bp++)=(fp++)->s;break;
			case 5: while(bp<be) *(bp++)=(fp++)->tau;break;
			case 6: while(bp<be) *(bp++)=(fp++)->chi;break;
			case 7: while(bp<be) *(bp++)=(fp++)->chi*TZ;break;
			case 8: while(bp<be) *(bp++)=(fp++)->dev();break;
			case 9: while(bp<be) *(bp++)=(fp++)->X;break;
			case 10: while(bp<be) *(bp++)=(fp++)->Y;break;
			case 11: while(bp<be) *(bp++)=(fp++)->cu;
		}
		if(!st) *bp=mode==9?buf[1]+bx-ax:buf[1];
		fwrite(buf,sizeof(float),l+1,outf);
	}

	// Close the file
	fclose(outf);
}

/** Computes the measure of the degree of quasi-staticity.
 * \param[out] Q the measure, integrated over the entire domain.
 * \param[out] max_qs the maximum value of the integrand.
 * \param[in] staggered whether to use the staggered velocity, only appropriate
 *			for quasi-static simulations. */
void shear_sim::qs_measure(double &Q,double &max_qs,bool staggered) {
	const double qsf1=6*K*mu*(K+(1/3.0)*mu),qsf2=6*K*mu*mu;
	Q=max_qs=0;
#pragma omp parallel for
	for(int j=1;j<n-1;j++) {
		double ux,uy,vx,vy,l=0,max_l=0,l1,l2,l3;

		// Evaluate the gradient of velocity
		for(int i=0;i<m;i++) {
			c_field *fp=fm+(i+m*j);
			if(staggered) {
				vec *vp=vel+(i+m*j);
				int d=i==m-1?1-m:1;
				ux=0.5*xsp*(-vp->x+vp[d].x-vp[m].x+vp[m+d].x),
				uy=0.5*ysp*(-vp->x-vp[d].x+vp[m].x+vp[m+d].x);
				vx=0.5*xsp*(-vp->y+vp[d].y-vp[m].y+vp[m+d].y);
				vy=0.5*ysp*(-vp->y-vp[d].y+vp[m].y+vp[m+d].y);
			} else {
				c_field &fr=*(fp+(i==m-1?1-m:1)),&fl=*(fp+(i==0?m-1:-1));
				c_field &fd=*(fp-m),&fu=*(fp+m);
				ux=0.5*xsp*(fr.u-fl.u);
				vx=0.5*xsp*(fr.v-fl.v);
				uy=0.5*ysp*(fu.u-fd.u);
				vy=0.5*ysp*(fu.v-fd.v);
			}
			l1=ux+vy;l2=ux-vy;l3=uy+vx;
			fp->cu=qsf1*l1*l1+qsf2*(l2*l2+l3*l3);
			l+=fp->cu;
			if(max_l<fp->cu) max_l=fp->cu;
		}

		// Store the computed values
#pragma omp critical
		{
			Q+=l;
			if(max_qs<max_l) max_qs=max_l;
		}
	}

	// Apply scaling to the measure
	Q=sqrt(Q*dx*dy);
}

/** Calculates the Hencky strain and stores it in the "cu" field entry. */
void shear_sim::compute_strain() {

	// Compute the stress on the regular grid
#pragma omp parallel for
	for(int j=0;j<n;j++) {
		c_field *fp=fm+j*m;
		for(int i=0;i<m;i++,fp++) {
			double l1,l2,J,Jinv,Xx,Xy,Yx,Yy;
			mat F,FT,FTF,sigma,Lam;

			// Calculate x derivatives of the reference map fields
			if(i==0) {
				Xx=(fp+1)->X-(fp+(m-1))->X+(bx-ax);
				Yx=(fp+1)->Y-(fp+(m-1))->Y;
			} else if(i==m-1) {
				Xx=(fp+(1-m))->X-(fp-1)->X-(bx-ax);
				Yx=(fp+(1-m))->Y-(fp-1)->Y;
			} else {
				Xx=(fp+1)->X-(fp-1)->X;
				Yx=(fp+1)->Y-(fp-1)->Y;
			}
			Xx*=0.5*xsp;Yx*=0.5*xsp;

			// Calculate the y derivatives of the reference map fields
			c_field *fu,*fd;
			double yfac;
			if(j==0) {yfac=ysp;fu=fp+m;fd=fp;}
			else if(j==n-1) {yfac=ysp;fu=fp+m;fd=fp;}
			else {yfac=0.5*ysp;fu=fp+m;fd=fp+m;}
			Xy=yfac*(fu->X-fd->X);
			Yy=yfac*(fu->Y-fd->Y);

			// Compute the deformation gradient tensor
			J=1/(Jinv=Xx*Yy-Xy*Yx);
			F=mat(Yy*J,-Xy*J,-Yx*J,Xx*J);

			// Compute the strain tensor
			FT=F.transpose();
			FTF=FT*F;
			FTF.sym_eigenvectors(l1,l2,Lam);
			if(l1>0||l2>0) {
				sigma=Lam.transpose()*mat(0.5*log(l1),0,0,0.5*log(l2))*Lam;
				fp->cu=sigma.devmod();
			} else fp->cu=0;
		}
	}
}

/** Initializes the effective temperature field as a bicubic interpolation of a
 * given array, scaling the array values.
 * \param[in] (mm,nn) the dimensions of the given array.
 * \param[in] arr a pointer to the array.
 * \param[in] lo an additive constant to scale by.
 * \param[in] sca an multiplicative constant to scale by. */
void shear_sim::initialize_chi_bicubic(int mm,int nn,double *arr,double lo,double sca) {
	int mmnn=mm*nn;
	double *xd=new double[3*mmnn],*yd=xd+mmnn,*xyd=yd+mmnn;

	// Calculate x derivative of the given array
#pragma omp parallel for
	for(int j=0;j<nn;j++) {
		int ij=j*mm;
		for(int i=0;i<mm;i++,ij++) xd[ij]=0.5*(arr[i==mm-1?ij-mm+1:ij+1]-arr[i==0?ij+mm-1:ij-1]);
	}

	// Calculate y and xy derivatives of the given array
#pragma omp parallel for
	for(int j=0;j<nn;j++) {
		int i,ij=j*mm;
		if(j==0) {
			for(i=0;i<mm;i++,ij++) {
				yd[ij]=arr[ij+mm]-arr[ij];
				xyd[ij]=xd[ij+mm]-xd[ij];
			}
		} else if (j==nn-1) {
			for(i=0;i<mm;i++,ij++) {
				yd[ij]=arr[ij]-arr[ij-mm];
				xyd[ij]=xd[ij]-xd[ij-mm];
			}
		} else {
			for(i=0;i<mm;i++,ij++) {
				yd[ij]=0.5*(arr[ij+mm]-arr[ij-mm]);
				xyd[ij]=0.5*(xd[ij+mm]-xd[ij-mm]);
			}
		}
	}

	// Assemble the effective temperature using bicubic interpolation
	double xsca=dx*mm/(bx-ax),ysca=dy*nn/(by-ay);
#pragma omp parallel for
	for(int j=0;j<n;j++) {
		double Y=(j+0.5)*ysca-0.5,mY;
		int jj=int(Y);
		if(jj<0) jj=0;else if(jj>nn-2) jj=nn-2;
		Y-=jj;mY=1-Y;
		double b0=Y*mY*mY,b1=mY*mY*(1+2*Y),b2=Y*Y*(3-2*Y),b3=-Y*Y*mY;
		for(int i=0;i<m;i++) {
			double X=(i+0.5)*xsca-0.5,mX;
			if(X<0) X+=mm;
			int ii=int(X),ij,ijp;
			if(ii>=mm) ii-=mm;
			ij=ii+mm*jj;ijp=(ii==mm-1?ii-mm:ii)+1+mm*jj;
			X-=ii;mX=1-X;
			double c0=X*mX*mX,c1=mX*mX*(1+2*X),c2=X*X*(3-2*X),c3=-X*X*mX;
			fm[i+m*j].chi=lo+sca*(b0*(xyd[ij]*c0+yd[ij]*c1+yd[ijp]*c2+xyd[ijp]*c3)
				+b1*(xd[ij]*c0+arr[ij]*c1+arr[ijp]*c2+xd[ijp]*c3)
				+b2*(xd[ij+mm]*c0+arr[ij+mm]*c1+arr[ijp+mm]*c2+xd[ijp+mm]*c3)
				+b3*(xyd[ij+mm]*c0+yd[ij+mm]*c1+yd[ijp+mm]*c2+xyd[ijp+mm]*c3));
		}
	}

	// Delete the temporary space
	delete [] xd;
}

/** Calculates the normalizing factors for the gaussian blurring.
 * \param[in] l the length scale of the gaussian, measured in gridpoints.
 * \param[in] cut the number of gridpoints to cut the PDF off at.
 * \return The normalization factor for the variance. */
double shear_sim::gaussian_normal_factor(double llinv,int cut) {
	double vfac=1,q,q2;
	for(int j=1;j<cut;j++) {
		q=exp(-j*j*llinv);
		for(int k=0;k<=j;k++) {
			q2=q*exp(-k*k*llinv);
			vfac+=(k==0||k==j?4:8)*q2*q2;
		}
	}
	printf("%g\n",vfac);
	return sqrt(vfac);
}

/** Initializes the simulation fields, setting the effective temperature field
 * to be a smoothed
 * \param[in] l the smoothing length scale, measured in gridpoints.
 * \param[in] tem_avg the average effective temperature.
 * \param[in] tem_stdev the standard deviation of the effective temperature. */
void shear_sim::initialize_random(double l,double tem_avg,double tem_stdev) {

	// Calculate cutoff and normalizing Gaussian factors
	int cut=int(l*gaussian_fac),nx=n+(cut<<1);
	double llinv=0.5/(l*l),chi0=tem_avg/TZ,chi1=tem_stdev/(gaussian_normal_factor(llinv,cut)*TZ);

	// Create grid of Gaussian noise
	double *no=new double[m*nx],*np=no,*ne=no+m*nx;
	for(;np<ne-1;np+=2) box_muller(*np,np[1]);
	if(np==ne-1) {
		double r1;
		box_muller(*np,r1);
	}

	// Initialize the simulation fields
	np=no+m*cut;
#pragma omp parallel for
	for(int j=0;j<=n;j++) {
		double Y=ay+j*dy;
		c_field *fp=fm+m*j;
		for(int i=0;i<m;i++,fp++) {
			c_field &f=*fp;
			double gau=0;

			// Set stresses and velocities to zero
			f.u=f.v=0;

			// Set reference map field to match initial coordinate
			f.X=ax+i*dx;f.Y=Y;

			// Calculate chi as a gaussian smoothing of the random field
			if(j<n) {
				for(int cj=j-cut;cj<=j+cut;cj++) for(int ci=i-cut;ci<=i+cut;ci++)
					gau+=np[cj*m+ci%m]*exp(-llinv*((cj-j)*(cj-j)+(ci-i)*(ci-i)));
				f.chi=chi0+gau*chi1;
				f.p=f.q=f.s=f.tau=0;
			}
		}
	}

	// Remove Gaussian noise grid
	delete [] no;
}

/** Calculates two normally distributed random numbers.
 * \param[out] (n0,n1) two normally distributed numbers. */
void shear_sim::box_muller(double &r0,double &r1) {
	double x,y,r;
	do {
		x=2.0*rand()/RAND_MAX-1;
		y=2.0*rand()/RAND_MAX-1;
		r=x*x+y*y;
	} while(r==0.0||r>1.0);
	double d=sqrt(-2.0*log(r)/r);
	r0=x*d;
	r1=y*d;
}

/** Calculates the tractions on the top and bottom surfaces.
 * \param[out] (bx,by) the components of the traction on the bottom surface.
 * \param[out] (tx,ty) the components of the traction on the top surface. */
void shear_sim::tractions(double &bx,double &by,double &tx,double &ty) {

	bx=by=tx=ty=0;

	// Calculate bottom traction, extrapolating the cell-centered stresses
	// to the boundary
	c_field *fp=fm,*fe=fp+m;
	while(fp<fe) {
		bx+=3*fp->tau-fp[m].tau;
		by+=3*fp->syy()-fp[m].syy();
		fp++;
	}

	// Calculate top traction
	fp=fm+(n-2)*m,fe=fp+m;
	while(fp<fe) {
		tx-=3*fp[m].tau-fp->tau;
		ty-=3*fp[m].syy()-fp->syy();
		fp++;
	}

	// Normalize by grid spacing, also taking into account factor of two in
	// extrapolation formula
	bx*=0.5*dx;by*=0.5*dx;
	tx*=0.5*dx;ty*=0.5*dx;
}
