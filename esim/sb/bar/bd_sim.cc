#include <cstring>

#include "common.hh"
#include "bd_sim.hh"
#include "extra_force.hh"
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
 * \param[in] chi_len_ the length scale for the chi diffusion.
 * \param[in] filter_stress_ the geometric stress filtering factor.
 * \param[in] t_scale_ a time scale associated with the plastic deformation.
 * \param[in] tmult_ a multiplier to apply to the default timestep size.
 * \param[in] adapt_fac_ a factor used in the adaptive timestepping of
 *			 plasticity.
 * \param[in] filename_ the filename of the output directory. */
bd_sim::bd_sim(const int m_,const int n_,const double ax_,const double bx_,
		const double ay_,const double by_,const double mu_,
		const double K_,const double viscosity_,const double chi_len_,
		const double filter_stress_,const double t_scale_,const double tmult_,
		const double adapt_fac_,
		stz_dynamics *stz_,const char *filename_,const unsigned int fflags_)
	: m(m_), n(n_), mn(m_*n_), me(m+1), ne(n+1), mne(me*ne), ax(ax_),
	ay(ay_), bx(bx_), by(by_), dx((bx_-ax_)/m_), dy((by_-ay_)/n_),
	xsp(1/dx), ysp(1/dy), mu(mu_), mu_inv(1/mu_), K(K_), viscosity(viscosity_),
	chi_len(chi_len_), filter_stress(filter_stress_), t_scale(t_scale_),
	tmult(tmult_), adapt_fac(adapt_fac_), stz(stz_), filename(filename_),
	fflags(fflags_), ls(m,n,ax+0.5*dx,ay+0.5*dy,dx,dy,-5*dx,5*dx),
	phi(ls.phi), c(ls.s), cc(ls.ss), fm((new c_field[mne])), vel(new vec[mne]),
	src(new vec[mne]), time(0), f_num(0), wall_bc(clamped), eforce(NULL),
	mgx(0), m_pad(22), qsm(NULL), mg(NULL), buf(new float[m>=63?m+2:64]), e_rm(fm),
	e_st(fm), e_reg(fm,m) {

	// Set default output bounds
	*obnd=0;obnd[1]=me;obnd[2]=0;obnd[3]=ne;
	obnd[4]=0;obnd[5]=m;obnd[6]=0;obnd[7]=n;

	// Set STZ related parameters
	TZ=stz->TZ;
	chi_inf=stz->chi_inf;

	// Set initial time for diagnostic purposes
	wc_time=wtime();
}

/** The class destructor frees the dynamically allocated memory. */
bd_sim::~bd_sim() {
	if(qsm!=NULL) {
		delete mg;
		delete qsm;
	}
	delete [] buf;
	delete [] src;
	delete [] vel;
	delete [] fm;
}

/** Initializes the simulation fields. */
void bd_sim::init_bar() {
	int j;
//	const double pi=3.1415926535897932384626433832795;
//	const double Bw=0.1,Bd=0.15;

	// Record initial wall position
	init_wallx=wallx;
	set_horizontal_bounds(dom_li,dom_ui);

	// Set up the level set to describe the object boundary
	double *phip=phi;
	for(j=0;j<n;j++) {
		double y=ay+(j+0.5)*dy,yy=fabs(y)-0.500001;
		for(int i=0;i<m;i++) {
			//x=ax+(i+0.5)*dx;
			//*(phip++)=fabs(x)<Bw?yy-Bd*(0.5+0.5*cos(x*pi/Bw)):yy;
			//*(phip++)=fabs(x)<Bw?max(yy,y-(0.500001-Bd*(0.5+0.5*cos(x*pi/Bw)))):yy; // Notched Initial Condition
			*(phip++)=yy; // Flat initial condition

		}
	}

	// Build the band, which will set the status array, and use it to set
	// the staggered status array
	ls.build_band();
	set_staggered_status();

	// Set the fields at points that are inside the object
	double chi0=620/TZ;
	double Dyy=strain_rate_yy();
#pragma omp parallel for
	for(j=0;j<=n;j++) {
		double y=ay+j*dy;
		for(int i=0;i<=m;i++) {
			int ij=i+m*j,ije=i+me*j;
			c_field &f=fm[ije];
			if(i<m&&j<n&&c[ij]<4) {
				f.p=f.q=f.s=f.tau=0;
				f.chi=chi0;
			}
			if(cc[ije]<=4) {
				double x=ax+i*dx;
				f.X=x;
				f.Y=y;
				f.v=Dyy*y;
				f.u=x>wallx?wallu:(x<-wallx?-wallu:wallu/wallx*x);
			}
		}
	}

	// Carry out initial field extrapolation
	extrapolate();
}

void bd_sim::set_staggered_status() {
#pragma omp parallel for
	for(int je=0;je<=n;je++) {
		int k,l,*cp=c+m*je,*ccp=cc+me*je;
		for(int ie=0;ie<=m;ie++,cp++,ccp++) {

			// Examine the four regular gridpoints adjacent to this
			// gridpoint, and count how many are in the fluid phase.
			k=l=0;
			if(ie>0) {
				if(je>0) {k++;if(cp[-m-1]>3) l++;}
				if(je<n) {k++;if(cp[-1]>3) l++;}
			}
			if(ie<m) {
				if(je>0) {k++;if(cp[-m]>3) l++;}
				if(je<n) {k++;if(*cp>3) l++;}
			}

			if(l==0) *ccp=0;
			else if(l==k) *ccp=7;
			else {
				double *phip=phi+ie+m*je,phit=0;
				if(ie>0) {
					if(je>0) phit+=phip[-m-1];
					if(je<n) phit+=phip[-1];
				}
				if(ie<m) {
					if(je>0) phit+=phip[-m];
					if(je<n) phit+=*phip;
				}
				*ccp=phit>0?7:0;
			}
		}
	}
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void bd_sim::solve(double duration,int frames) {
	double t_start=time,time_interval=duration/frames,target_time,t0,t1,t2;
	int l;
	const double dt=dx*dx*tmult;

	// Output the initial fields and record initial time
	if(f_num==0) {write_files(0);puts("# Output frame 0");}
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
		write_files(k+f_num);

		// Print diagnostic information
		t2=wtime();
		printf("# Output frame %d, t=%g [%d, %.8g s, %.8g s]\n",k,time,l,t1-t0,t2-t1);
		t0=t2;
	}
	f_num+=frames;
}

/** Carries out the simulation for a specified time interval using the
 * quasi-static simulation method, periodically saving the output.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save.
 * \param[in] steps the number of quasi-static steps to take per frame. */
void bd_sim::solve_quasistatic(double duration,int frames,int steps) {
	double time_interval=duration/frames,t0,t1,t2;
	double isteps=1.0/steps;

	// Output the initial fields
	if(f_num==0) {
		write_files(0);
		puts("# Output frame 0");
	}
	t0=wtime();

	for(int k=1;k<=frames;k++) {

		for(int j=0;j<steps;j++) step_forward_quasistatic(time_interval*isteps);
		t1=wtime();

		// Output the fields
		write_files(k+f_num);

		// Print diagnositic information
		t2=wtime();
		printf("# Output frame %d t=%g [%d, %.8g s, %.8g s]\n",k,time,steps,t1-t0,t2-t1);
		t0=t2;
	}
	f_num=frames;
}

/** Prototype adaptive quasi-static solver. */
void bd_sim::solve_quasistatic_adaptive(double duration,int frames,int steps) {
	double time_interval=duration/frames,t0,t1,t2;
	double isteps=1.0/steps;

	// Output the initial fields
	if(f_num==0) {
		write_files(0);
		puts("# Output frame 0");
	}
	t0=wtime();

	for(int k=1;k<=frames;k++) {

		for(int j=0;j<steps;j++) step_forward_quasistatic(time_interval*isteps);
		t1=wtime();

		// Output the fields
		write_files(k+f_num);

		// Print diagnositic information
		t2=wtime();
		printf("# Output frame %d t=%g [%d, %.8g s, %.8g s]\n",k,time,steps,t1-t0,t2-t1);
		t0=t2;
	}
	f_num=frames;
}

/** Steps the simulation fields forward using the direct update procedure.
 * \param[in] dt the time step to use. */
void bd_sim::step_forward(double dt) {
	time+=dt;

	// Wall setup
	wall_setup(dt);

	// Direct step and post-processing
	direct_step(dt);post_process();
}

/** Steps the simulation fields forward using the quasistatic update procedure.
 * \param[in] dt the time step to use. */
void bd_sim::step_forward_quasistatic(double dt) {

	// Update the simulation time
	time+=dt;
	if(eforce!=NULL) eforce->compute(dt);

	// Wall setup
	wall_setup(dt);set_boundaries();

	// Carry out the intermediate step that considers the advective and
	// plasiticty terms
	advection_step(dt);post_process(false);

	// Apply the projection step to ensure quasistaticity
	projection_step(dt);post_process();

	if(eforce!=NULL) eforce->update();
}

void bd_sim::wall_setup(double dt) {
	wallx+=wallu*dt;
	int l,o;

	set_horizontal_bounds(l,o);
	if(l<dom_li) simple_extrapolation(l,l+1,l+2);
	dom_li=l;

	if(o>dom_ui) simple_extrapolation(o-1,o-2,o-3);
	dom_ui=o;
}

void bd_sim::simple_extrapolation(int i,int i2,int i3) {
	int jm;
	c_field *fp=fm+i,*fp2=fm+i2,*fp3=fm+i3;
	for(jm=0;jm<mne;jm+=me) {
		fp[jm].p=2*fp2[jm].p-fp3[jm].p;
		fp[jm].q=2*fp2[jm].q-fp3[jm].q;
		fp[jm].s=2*fp2[jm].s-fp3[jm].s;
		fp[jm].tau=2*fp2[jm].tau-fp3[jm].tau;
		fp[jm].chi=2*fp2[jm].chi-fp3[jm].chi;
	}
}

void bd_sim::post_process(bool all_fields) {
	if(dom_li>0) simple_extrapolation(dom_li-1,dom_li,dom_li+1);
	if(dom_ui<m) simple_extrapolation(dom_ui,dom_ui-1,dom_ui-2);
	extrapolate(all_fields);
	set_boundaries();
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void bd_sim::direct_step(double dt) {
	const double third=1./3;
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt;
	double hxxv=2*xsp*hx*viscosity,hyyv=2*ysp*hy*viscosity,hcv=2*(hxxv+hyyv);

	// Move the level set
	ls.move(this,dt);
	set_staggered_status();

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		c_field fg[4];
		for(int i=0,*cp=c+j*m,*ccp=cc+j*me;i<m;i++,cp++,ccp++) {
			c_field *fp=fm+i+me*j,&f=*fp;

			// Calculate the update to the reference map and
			// velocity
			if(j>0&&i>0&&*ccp<4) {
				double ux,vx,Xx,Yx,fx,uy,vy,Xy,Yy,fy,uc=f.u,vc=f.v,uvisc,vvisc;
				c_field &fd=fp[-me],&fl=fp[-1],&fr=fp[1],&fu=fp[me];

				// Calculate ENO derivatives of regular fields
				uc>0?(i>1&&ccp[-2]<4?rmv_eno2(ux,vx,Xx,Yx,hx,fr,f,fl,fp[-2])
						    :rmv_one_sided(ux,vx,Xx,Yx,hx,f,fl))
				    :(i<m-1&&ccp[2]<4?rmv_eno2(ux,vx,Xx,Yx,-hx,fl,f,fr,fp[2])
						     :rmv_one_sided(ux,vx,Xx,Yx,-hx,f,fr));
				vc>0?(j>1&&ccp[-2*m]<4?rmv_eno2(uy,vy,Xy,Yy,hy,fu,f,fd,fp[-2*me])
						      :rmv_one_sided(uy,vy,Xy,Yy,hy,f,fd))
				    :(j<n-1&&ccp[2*m]<4?rmv_eno2(uy,vy,Xy,Yy,-hy,fd,f,fu,fp[2*me])
						       :rmv_one_sided(uy,vy,Xy,Yy,-hy,f,fu));

				// Calculate net force due to stress imbalance
				net_force(fp,i,j,fx,fy,-hx,-hy,fg);

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

			if(*cp<4) {

				// Calculate deformation rate, staggered
				// velocity, and spin
				double uc=0.25*(f.u+fp[1].u+fp[me].u+fp[me+1].u),
				       vc=0.25*(f.v+fp[1].v+fp[me].v+fp[me+1].v),

				// Calculate the deformation rate
				       ux=hx*(-f.u+fp[1].u-fp[me].u+fp[me+1].u),
				       uy=hy*(-f.u-fp[1].u+fp[me].u+fp[me+1].u),
				       vx=hx*(-f.v+fp[1].v-fp[me].v+fp[me+1].v),
				       vy=hy*(-f.v-fp[1].v+fp[me].v+fp[me+1].v),
				       two_omega=vx-uy,dtmp=f.dev();

				// Calculate ENO derivatives of staggered fields
				st_field f1,f2;bool dd,uu;
				c_field &fl=(dd=bc_deriv(i,j,i-1,j,*fg))?fp[-1]:*fg,
					&fr=(uu=bc_deriv(i,j,i+1,j,fg[1]))?fp[1]:fg[1];
				f1=uc>0?(dd&&i>1&&cp[-2]<4?st_eno2(hx,fr,f,fl,fp[-2]):st_one_sided(hx,f,fl))
				       :(uu&&i<m-2&&cp[2]<4?st_eno2(-hx,fl,f,fr,fp[2]):st_one_sided(-hx,f,fr));
				c_field	&fd=(dd=bc_deriv(i,j,i,j-1,*fg))?fp[-me]:*fg,
					&fu=(uu=bc_deriv(i,j,i,j+1,fg[1]))?fp[me]:fg[1];
				f2=vc>0?(dd&&j>1&&cp[-2*m]<4?st_eno2(hy,fu,f,fd,fp[-2*me]):st_one_sided(hy,f,fd))
				       :(uu&&j<n-2&&cp[2*m]<4?st_eno2(-hy,fd,f,fu,fp[2*me]):st_one_sided(-hy,f,fu));

				// Calculate the adaptive effective temperature
				// update
				f.cchi=f.chi;
				dtmp=adaptive_plastic_term(dtmp,f.cchi,f.ddev,dt);
				f.cchi-=uc*f1.chi+vc*f2.chi;

				// Calculate updates to stress
				f.cp=-uc*f1.p-vc*f2.p-K*(ux+vy);
				f.cq=-uc*f1.q-vc*f2.q-mu*third*(ux+vy)-dtmp*f.q;
				f.cs=-uc*f1.s-vc*f2.s-two_omega*f.tau+mu*(ux-vy)-dtmp*f.s;
				f.ctau=-uc*f1.tau-vc*f2.tau+two_omega*f.s+mu*(vx+uy)-dtmp*f.tau;
			}
		}
	}

	// Add in the chi diffusion
	chi_diffusion(dt);

	// Apply the updates to the fields
#pragma omp parallel for
	for(int j=0;j<n;j++) {
		int ij=m*j,ijf=ij+m;
		if(j>0) {
			if(c[ij]<4) fm[ij+j].update_staggered();
			ij++;
			for(;ij<ijf;ij++) {
				if(c[ij]<4) fm[ij+j].update_staggered();
				if(cc[ij+j]<4) fm[ij+j].update_regular();
			}
		} else for(;ij<ijf;ij++) if(c[ij]<4) fm[ij+j].update_staggered();
	}
}

/** Carries out an intermediate step forward, by evaluating the advective and
 * plasticity terms in the equations only.
 * \param[in] dt the time step to use. */
void bd_sim::advection_step(double dt) {
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt;
//	double hxxv=2*xsp*hx,hyyv=2*ysp*hy,
	double hxxv=filter_stress,hyyv=filter_stress,hcv=2*(hxxv+hyyv);
	dom_lj=n;dom_uj=0;

	// Move the level set
	ls.move(this,dt);
	set_staggered_status();

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		bool in_row=false;
		c_field fg[4];
		for(int i=dom_li,*cp=c+j*m+dom_li,*ccp=cc+j*me+dom_li;i<dom_ui;i++,cp++,ccp++) {
			c_field *fp=fm+i+me*j,&f=*fp;

			// Calculate the update to the reference map
			if(j>0&&i>dom_li&&*ccp<4) {
				double Xx,Yx,Xy,Yy,uc=f.u,vc=f.v;
				c_field &fd=fp[-me],&fl=fp[-1],&fr=fp[1],&fu=fp[me];

				// Calculate ENO derivatives of regular fields
				uc>0?(i>dom_li+1&&ccp[-2]<4?rm_eno2(Xx,Yx,hx,fr,f,fl,fp[-2])
						    :rm_one_sided(Xx,Yx,hx,f,fl))
				    :(i<dom_ui-1&&ccp[2]<4?rm_eno2(Xx,Yx,-hx,fl,f,fr,fp[2])
						     :rm_one_sided(Xx,Yx,-hx,f,fr));
				vc>0?(j>1&&ccp[-2*m]<4?rm_eno2(Xy,Yy,hy,fu,f,fd,fp[-2*me])
						      :rm_one_sided(Xy,Yy,hy,f,fd))
				    :(j<n-1&&ccp[2*m]<4?rm_eno2(Xy,Yy,-hy,fd,f,fu,fp[2*me])
						       :rm_one_sided(Xy,Yy,-hy,f,fu));

				// Calculate update to reference map fields
				f.cX=-uc*Xx-vc*Xy;
				f.cY=-uc*Yx-vc*Yy;
			}

			if(*cp<4) {
				in_row=true;

				// Calculate deformation rate, staggered
				// velocity, and spin
				double uc=0.25*(f.u+fp[1].u+fp[me].u+fp[me+1].u),
				       vc=0.25*(f.v+fp[1].v+fp[me].v+fp[me+1].v),
				       two_omega=hx*(-f.v+fp[1].v-fp[me].v+fp[me+1].v)
					        +hy*(-f.u-fp[1].u+fp[me].u+fp[me+1].u);

				// Calculate ENO derivatives of staggered fields
				st_field f1,f2;bool dd,uu;
				c_field &fl=(dd=bc_deriv_proj(i,j,i-1,j,*fg))?fp[-1]:*fg,
					&fr=(uu=bc_deriv_proj(i,j,i+1,j,fg[1]))?fp[1]:fg[1];
				f1=uc>0?(dd&&i>1&&cp[-2]<4?st_eno2(hx,fr,f,fl,fp[-2]):st_one_sided(hx,f,fl))
				       :(uu&&i<m-2&&cp[2]<4?st_eno2(-hx,fl,f,fr,fp[2]):st_one_sided(-hx,f,fr));
				c_field	&fd=(dd=bc_deriv_proj(i,j,i,j-1,*fg))?fp[-me]:*fg,
					&fu=(uu=bc_deriv_proj(i,j,i,j+1,fg[1]))?fp[me]:fg[1];
				f2=vc>0?(dd&&j>1&&cp[-2*m]<4?st_eno2(hy,fu,f,fd,fp[-2*me]):st_one_sided(hy,f,fd))
				       :(uu&&j<n-2&&cp[2*m]<4?st_eno2(-hy,fd,f,fu,fp[2*me]):st_one_sided(-hy,f,fu));

				// Calculate filtering terms
				double pfil=hxxv*(fr.p+fl.p)+hyyv*(fu.p+fd.p)-hcv*f.p,
				       qfil=hxxv*(fr.q+fl.q)+hyyv*(fu.q+fd.q)-hcv*f.q,
				       sfil=hxxv*(fr.s+fl.s)+hyyv*(fu.s+fd.s)-hcv*f.s,
				       taufil=hxxv*(fr.tau+fl.tau)+hyyv*(fu.tau+fd.tau)-hcv*f.tau,dtmp=f.dev();

				// Calculate the adaptive effective temperature
				// update
				f.cchi=f.chi;
				dtmp=adaptive_plastic_term(dtmp,f.cchi,f.ddev,dt);
				f.cchi-=uc*f1.chi+vc*f2.chi;

				// Calculate updates to stress
				f.cp=-uc*f1.p-vc*f2.p+pfil;
				f.cq=-uc*f1.q-vc*f2.q-dtmp*f.q+qfil;
				f.cs=-uc*f1.s-vc*f2.s-two_omega*f.tau-dtmp*f.s+sfil;
				f.ctau=-uc*f1.tau-vc*f2.tau+two_omega*f.s-dtmp*f.tau+taufil;
			}
		}

		// Set vertical bounds
		if(in_row) {
#pragma omp critical
			{
				if(j<dom_lj) dom_lj=j;
				if(j>dom_uj) dom_uj=j;
			}
		}
	}

	// Check for errors
	if(dom_lj==n) {
		fputs("No gridpoints\n",stderr);
		exit(1);
	}

	// Adjust counters to account for ghost layers needed in projection
	// step
	dom_lj-=1;dom_uj+=2;

	// Add in the chi diffusion
	chi_diffusion(dt);

	// Apply the updates to the fields
#pragma omp parallel for
	for(int j=dom_lj;j<dom_uj;j++) {
		int ij=m*j+dom_li,ijf=m*j+dom_ui;
		if(j>dom_li) {
			if(c[ij]<4) fm[ij+j].update_staggered();
			ij++;
			for(;ij<ijf;ij++) {
				if(c[ij]<4) fm[ij+j].update_staggered();
				if(cc[ij+j]<4) fm[ij+j].update_ref_map();
			}
		} else for(;ij<ijf;ij++) if(c[ij]<4) fm[ij+j].update_staggered();
	}
}

/** Adds in contributions to chi diffusion to change in the chi field. */
void bd_sim::chi_diffusion(double dt) {

	// Calculate normalizing factor to apply when computing the diffusion.
	// The stored ddev values contain dt*2*mu*Dpl. The 2*mu must be scaled
	// out. There will be an additional factor of 2 because Dpl averages on
	// edges are computed from two numbers. Grid spacings of dx*dx/dt must
	// also be dealt with.
	const double dfac=0.25*chi_len*chi_len*mu_inv,dfacx=dfac*xsp*xsp,dfacy=dfac*ysp*ysp;

#pragma omp parallel for
	for(int j=0;j<n;j++) for(int i=0;i<m;i++) {
		int *cp=c+i+m*j;
		if(*cp<4) {

			c_field *fp=fm+(i+me*j),&f=*fp,&fl=fp[-1],&fr=fp[1],&fd=fp[-me],&fu=fp[me];
			if(cp[-1]<4&&cp[1]<4) {

				// Compute D_chi terms
				double Dl=f.ddev+fl.ddev,Dr=f.ddev+fr.ddev;

				// Diffusion mediated by Dpl
				f.cchi+=dfacx*(Dr*fr.chi-(Dr+Dl)*f.chi+Dl*fl.chi);
			}
			if(cp[-m]<4&&cp[m]<4) {

				// Compute D_chi terms
				double Dd=f.ddev+fd.ddev,Du=f.ddev+fu.ddev;

				// Diffusion mediated by Dpl
				f.cchi+=dfacy*(Du*fu.chi-(Du+Dd)*f.chi+Dd*fd.chi);
			}
		}
	}
}

/** Carries out a projection step of the stresses using a multigrid solve,
 * to enforce the quasistaticity constraint.
 * \param[in] dt the timestep to use. */
void bd_sim::projection_step(double dt) {
	const double third=1/3.0;
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt;

	// If the bounds have changed from last time, then delete the previous
	// multigrid setup class
	if(qsm!=NULL&&(dom_li!=qsm_li||dom_ui!=qsm_ui||dom_lj!=qsm_lj||dom_uj!=qsm_uj)) {
		printf("New multigrid setup: (%d,%d) (%d,%d) [%d]\n",dom_li,dom_ui,dom_lj,dom_uj,m_pad);
		delete qsm;qsm=NULL;
		delete mg;
	}

	// Allocate the multigrid setup class if needed
	if(qsm==NULL) {
		qsm_li=dom_li;qsm_ui=dom_ui;
		qsm_lj=dom_lj;qsm_uj=dom_uj;
		qsm=new qs_multi_bd(*this,qsm_li,qsm_ui,qsm_lj-m_pad,qsm_uj+m_pad);
		mg=new tgmg<qs_multi_bd,vec,mat>(*qsm,src,vel);
		mg->verbose=1;
	}

	// Adjust multigrid pad if previous time was not effective
	else if(mg->conv_rate<5e-2) {
		adjust_multigrid_pad();
		printf("Adjust padding to %d to try and get better convergence\n",m_pad);
	}

	// Initialize constants in the algebraic system
	qsm->init(viscosity,dt);
	//qsm->adjust_ro(1e-8*pow(1e6,0.01*mgx));

	// Set up the algebraic systems in the grid hierarchy. Reduce the
	// tuning factor if needed.
	mg->setup();
	if(mgx>0) mgx--;

	// First try a few solves, and adjust the padding if they don't work
	int count=0;
	while(count<5&&!solve_linear_system()) {
		adjust_multigrid_pad();
		qsm->init(viscosity,dt);
		mg->setup();
		mg->clear_z();
		printf("Multigrid failed; adjust padding to %d\n",m_pad);
		count++;
	}

	// The previous attempts weren't able to solve the system. Now try
	// lowering the number of multigrid levels
	if(count==5) {
		while(!solve_linear_system()) {
			mg->clear_z();
			mgx+=12;
			//qsm->adjust_ro(1e-8*pow(1e6,0.01*mgx));
			//mg->setup();
			printf("Multigrid failed; increase tuning factor to %d\n",mgx);
		}
	}

	// Apply the correction to the stress tensor, based upon the computed
	// velocity. In addition, set up a representation of the velocity on
	// the primary grid
	int qm=qsm->m;
#pragma omp parallel for
	for(int j=dom_lj;j<dom_uj;j++) for(int i=dom_li;i<dom_ui;i++) {
		int ij=i+m*j,ije=ij+j;
		c_field &f=fm[ije];
		vec *vp=vel+(j-dom_lj+m_pad)*qm+(i-dom_li);
		if(j>0&&i>0) {//&&cc[ije]<4) {
			if(cc[ije]<400) {f.u=vp->x;f.v=vp->y;}
			else f.u=f.v=0;
		}
		if(c[ij]<4) {
			double ux=hx*(-vp->x+vp[1].x-vp[qm].x+vp[qm+1].x),
			       uy=hy*(-vp->x-vp[1].x+vp[qm].x+vp[qm+1].x),
			       vx=hx*(-vp->y+vp[1].y-vp[qm].y+vp[qm+1].y),
			       vy=hy*(-vp->y-vp[1].y+vp[qm].y+vp[qm+1].y);
			f.p+=-K*(ux+vy);
			f.q+=-mu*third*(ux+vy);
			f.s+=mu*(ux-vy);
			f.tau+=mu*(uy+vx);
		}
	}
}

/** Computes ghost values for the fields at the boundaries.
 * \param[out] (uu,vu,pu,qu,su,tauu,chiu) the computed ghost values.
 * \param[in] (i,j) the current gridpoint.
 * \param[in] (ni,nj) the gridpoint at which to compute the ghost values
 * \return True if this is a regular gridpoint, false otherwise. */
bool bd_sim::bc_deriv(int i,int j,int ni,int nj,c_field &fg) {
	if (ni<0||ni>=m||nj<0||nj>=n) {
		c_field &f=fm[i+me*j],&fb=fm[2*i-ni+me*(2*j-nj)];
		fg.p=2*f.p-fb.p;fg.q=2*f.q-fb.q;
		fg.s=2*f.s-fb.s;fg.tau=2*f.tau-fb.tau;
		fg.chi=2*f.chi-fb.chi;
		return false;
	}
	int nij;
	if(c[nij=ni+m*nj]<4) return true;
	int ij=i+m*j;
	c_field &f=fm[ij+j],&fn=fm[nij+nj];
	double z=phi[ij]/(phi[ij]-phi[nij]),
	       pp=f.p*(1-z)+z*fn.p,qp=f.q*(1-z)+z*fn.q,
	       sp=f.s*(1-z)+z*fn.s,taup=f.tau*(1-z)+z*fn.tau,
	       phix=phi_x(i,ij)*(1-z)+z*phi_x(ni,nij),
	       phiy=phi_y(ij)*(1-z)+z*phi_y(nij);
	project_stress(pp,qp,sp,taup,phix,phiy);
	int bi=2*i-ni,bj=2*j-nj;
	bool bin=bi>=0&&bi<m&&bj>=0&&bj<n;
	c_field &fb=bin?fm[bi+me*bj]:f;
	z=bin?2/(1+z):1/z;
	fg.p=fb.p*(1-z)+z*pp;
	fg.q=fb.q*(1-z)+z*qp;
	fg.s=fb.s*(1-z)+z*sp;
	fg.tau=fb.tau*(1-z)+z*taup;
	fg.chi=f.chi;
	return false;
}

/** Computes ghost values for the fields at the boundaries.
 * \param[out] (uu,vu,pu,qu,su,tauu,chiu) the computed ghost values.
 * \param[in] (i,j) the current gridpoint.
 * \param[in] (ni,nj) the gridpoint at which to compute the ghost values
 * \return True if this is a regular gridpoint, false otherwise. */
bool bd_sim::bc_deriv_proj(int i,int j,int ni,int nj,c_field &fg) {
	if (ni<dom_li||ni>=dom_ui||nj<0||nj>=n) {
		c_field &f=fm[i+me*j],&fb=fm[2*i-ni+me*(2*j-nj)];
		fg.p=2*f.p-fb.p;fg.q=2*f.q-fb.q;
		fg.s=2*f.s-fb.s;fg.tau=2*f.tau-fb.tau;
		fg.chi=2*f.chi-fb.chi;
		return false;
	} else return true;
}

/** Projects the components of stress to only have a sigma_tt component. If the
 * unnormalized derivative of the level set function is too close to zero, no
 * projection is carried out.
 * \param[in,out] (pp,qp,sp,taup) the components of the stress tensor.
 * \param[in] (phix,phiy) the unnormalized derivative of the level set function. */
void bd_sim::project_stress(double &pp,double &qp,double &sp,double &taup,double phix,double phiy) {
	const double third=1/3.0;
	double magn=phix*phix+phiy*phiy;
	if(magn<1e-16) return;
	magn=1/sqrt(magn);
	phix*=magn;phiy*=magn;
	double p2p=pp+qp,s33p=-pp+2*qp,
	       sin2t=2*phix*phiy,cos2t=1-2*phiy*phiy,
	       sigmap=-p2p-sp*cos2t-taup*sin2t;
	p2p=-0.5*sigmap;
	sp=-0.5*sigmap*cos2t;
	taup=-0.5*sigmap*sin2t;
	pp=third*(2*p2p-s33p);
	qp=third*(p2p+s33p);
}

/** Calculates the net force at a regular grid point using the stress tensors
 * on the four neighboring points of the staggered grids.
 * \param[in] fp a pointer to the grid point to consider.
 * \param[in] ij the regular grid point index.
 * \param[out] (fx,fy) the components of force.
 * \param[in] (xf,yf) prefactors to apply to the x and y derivatives that go into the computation of force.
 * \param[in] fg a pointer to ghost value storage. */
void bd_sim::net_force(c_field *fp,int i,int j,double &fx,double &fy,double xf,double yf,c_field* fg) {
	int ij=i+m*j;
	double *phip=phi+ij,phic=0.25*(*phip+phip[-1]+phip[-m]+phip[-m-1]);

	// Look at each corner, and construct ghost values if needed
	c_field &f0=diag_line(fp-me-1,ij-m-1,i>1&&j>1,phic,1,m,me,*fg),
		&f1=diag_line(fp-me,ij-m,i<m-1&&j>1,phic,-1,m,me,fg[1]),
		&f2=diag_line(fp-1,ij-1,i>1&&j<n-1,phic,1,-m,-me,fg[2]),
		&f3=diag_line(fp,ij,i<m-1&&j<n-1,phic,-1,-m,-me,fg[3]);

	// Compute the net force, applying the given prefactors
	fx=(-f0.p+f1.p-f2.p+f3.p-f0.q+f1.q-f2.q+f3.q+f0.s-f1.s+f2.s-f3.s)*xf
	  +(f0.tau+f1.tau-f2.tau-f3.tau)*yf;
	fy=(f0.tau-f1.tau+f2.tau-f3.tau)*xf
	  +(-f0.p-f1.p+f2.p+f3.p-f0.q-f1.q+f2.q+f3.q-f0.s-f1.s+f2.s+f3.s)*yf;
}

inline c_field& bd_sim::diag_line(c_field *fp,int ij,bool look_more,double phic,int d1,int d2,int d2e,c_field &fg) {
	double z,fz,phix,phiy,*phip=phi+ij;

	// If the corner is inside, then just make use of those values
	if(c[ij]<4) {
		if(!look_more||c[ij-d1-d2]<4) return *fp;
		z=*phip/(*phip-phip[-d1-d2]);
		if(z>0.5) return *fp;
		fz=1-z;
		phix=fz*(phip[-d1-d2]-phip[-d2])+z*(phip[-d1]-*phip);
		phiy=fz*(phip[-d1-d2]-phip[-d1])+z*(phip[-d2]-*phip);
		z+=1;fz=1-z;
	} else {
		z=(*phip*0.5-phic)/(*phip-phic);
		fz=1-z;
		phix=fz*(phip[d2]-phip[d1+d2])+z*(*phip-phip[d1]);
		phiy=fz*(phip[d1]-phip[d1+d2])+z*(*phip-phip[d2]);
	}

//	double xpos=(ij%m)+(d1==-1?-0.5:0.5),ypos=(ij/m)+(d2==-m?-0.5:0.5),gz=z-0.5;
	phix*=d1>0?-xsp:xsp;phiy*=d2>0?-ysp:ysp;
//	printf("%g %g %g %g %g %g\n",xpos,ypos,(d1==-1?gz:-gz),(d2==-m?gz:-gz),phix,phiy);

	// Compute the level set position along the diagonal. Evaluate
	// the normal vector and stress tensor at this position.
	double pp=z*fp->p+fz*fp[d1+d2e].p,qp=z*fp->q+fz*fp[d1+d2e].q,
	       sp=z*fp->s+fz*fp[d1+d2e].s,taup=z*fp->tau+fz*fp[d1+d2e].tau;

	// Project out normal components of stress
	project_stress(pp,qp,sp,taup,phix,phiy);

	// Calculate ghost stress based on extrapolation of the projected
	// stress to the corner
	z=1/z;
	c_field &fb=fp[d1+d2e];
	fg.p=fb.p*(1-z)+pp*z;
	fg.q=fb.q*(1-z)+qp*z;
	fg.s=fb.s*(1-z)+sp*z;
	fg.tau=fb.tau*(1-z)+taup*z;
	return fg;
}

/** Adaptively calculates the effect of plasticity on the shear stress
 * components and the effective temperature, by dividing up the timestep
 * into smaller intervals to ensure that the change in shear stress remains
 * small.
 * \param[in] sbar the magnitude of the shear stress at a grid point.
 * \param[in,out] chiv the effective temperature at the grid point.
 * \param[out] ddev a scaled measure of Dpl, needed for chi diffusion.
 * \param[in] dt the timetep.
 * \return The factor to scale the shear stresses by. */
double bd_sim::adaptive_plastic_term(double sbar,double &chiv,double &ddev,double dt) {
	const double small_dev_cutoff=1e-12;
	const double theta=0;
	const bool no_chidot=false;
	bool adapt=true;
	double dtl=dt,osbar=sbar,dplas,dchi1,dchi2,adt;

	// To avoid division-by-zero errors, check for the case when the
	// deviatoric stress is small, and just return with no changes
	if(sbar<small_dev_cutoff) {ddev=0;return 0;}
	do {

		// Calculate the plastic deformation
		dplas=2*t_scale*mu*stz->Dplastic(sbar,chiv,dchi1,dchi2);

		// If the plastic deformation is large, then adaptively
		// choose the timestep so as not to make a large alteration in
		// the deviatoric stress magnitude
		if(fabs(dtl*dplas)>adapt_fac) {
			adt=adapt_fac/fabs(dplas);
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
	ddev=osbar-sbar;
	return ddev/osbar;
}

/** Sets the fields in the top and bottom boundaries according to the boundary
 * conditions. */
void bd_sim::set_boundaries() {

	// Set bottom boundary
	for(int i=0;i<=m;i++) if(cc[i]<8) {
		fm[i].u=0;
		fm[i].v=-0.000;
		fm[i].X=ax+i*dx;
		fm[i].Y=ay;
	}

	// Set left and right boundaries
	double Dyy=strain_rate_yy();
	for(int j=1;j<n;j++) {
		int jme=j*me;
		double x,y=ay+j*dy;

		// Set left BCs
		for(int i=0;i<=dom_li;i++) {
			x=ax+i*dx;
			fm[jme+i].u=-wallu;
			fm[jme+i].v=y*Dyy;
			fm[jme+i].X=x+(wallx-init_wallx);
			fm[jme+i].Y=y;
		}

		// Set right BCs
		for(int i=dom_ui;i<=m;i++) {
			x=ax+i*dx;
			fm[jme+i].u=wallu;
			fm[jme+i].v=y*Dyy;
			fm[jme+i].X=x-(wallx-init_wallx);
			fm[jme+i].Y=y;
		}
	}

	// Set top boundary
	for(int i=0,ij=me*n;i<=m;i++,ij++) if(cc[ij]<8) {
		fm[ij].u=0;
		fm[ij].v=0;
		fm[ij].X=ax+i*dx;
		fm[ij].Y=by;
	}
}

/** Computes the yy-component of strain rate in the bar, assuming an elastic
 * response */
double bd_sim::strain_rate_yy() {

	if(wall_bc==clamped) {

		// Clamped boundary condition - there is no vertical motion
		return 0;
	} else if (wall_bc==co_thin) {

		// Scale factor, equivalent to nu/(1-nu)
		const double factor=(3*K-2*mu)/(3*K+4*mu);

		// Compute strain rate component. I think there are other things that
		// might be reasonable here.
		return -factor*wallu/wallx;
	} else {

		// Calculate average du/dx over the two end strips
		int o=0;
		double ux=0;
		for(int j=0;j<n;j++) {
		//for(int j=n/2-1;j<n/2+2;j++) { // Average just over a few gridpoints in the center of the bar
			int *ccp=cc+j*me;
			c_field *fp=fm+me*j;
			if(ccp[dom_li+4]==0&&ccp[dom_li+8]==0) {
				ux+=fp[dom_li+8].u-fp[dom_li+4].u;
				o++;
			}
			if(ccp[dom_ui-8]==0&&ccp[dom_ui-4]==0) {
				ux+=fp[dom_ui-4].u-fp[dom_ui-8].u;
				o++;
			}
		}
		if(o==0) {
			fputs("No gridpoints to compute du/dx\n",stderr);
			exit(1);
		}
		return -xsp*ux/o*0.25;
	}
}

/** Calculates one-sided derivatives of the fields using first-order
 * upwinding.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f1,f2) the fields to compute the derivative with.
 * \return The computed derivative. */
inline void bd_sim::rmv_one_sided(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f1,c_field &f2) {
	ud=2*hs*(f1.u-f2.u);
	vd=2*hs*(f1.v-f2.v);
	rm_one_sided(Xd,Yd,hs,f1,f2);
}

/** Calculates one-sided derivatives of the fields using first-order
 * upwinding.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f1,f2) the fields to compute the derivative with.
 * \return The computed derivative. */
inline void bd_sim::rm_one_sided(double &Xd,double &Yd,double hs,c_field &f1,c_field &f2) {
	hs*=2;
	Xd=hs*(f1.X-f2.X);
	Yd=hs*(f1.Y-f2.Y);
}

/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \return The computed derivative. */
inline void bd_sim::rmv_eno2(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3) {
	ud=hs*eno2(f0.u,f1.u,f2.u,f3.u);
	vd=hs*eno2(f0.v,f1.v,f2.v,f3.v);
	rm_eno2(Xd,Yd,hs,f0,f1,f2,f3);
}

/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \return The computed derivative. */
inline void bd_sim::rm_eno2(double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3) {
	Xd=hs*eno2(f0.X,f1.X,f2.X,f3.X);
	Yd=hs*eno2(f0.Y,f1.Y,f2.Y,f3.Y);
}

/** Calculates one-sided derivatives of the fields using first-order
 * upwinding.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f1,f2) the fields to compute the derivative with.
 * \return The computed derivative. */
inline st_field bd_sim::st_one_sided(double hs,c_field &f1,c_field &f2) {
	hs*=2;
	return st_field(hs*(f1.p-f2.p),hs*(f1.q-f2.q),hs*(f1.s-f2.s),hs*(f1.tau-f2.tau),
		hs*(f1.chi-f2.chi));
}

/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme, applying the shift to the X terms.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \return The computed derivative. */
st_field bd_sim::st_eno2(double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3) {
	return st_field(hs*eno2(f0.p,f1.p,f2.p,f3.p),hs*eno2(f0.q,f1.q,f2.q,f3.q),
		hs*eno2(f0.s,f1.s,f2.s,f3.s),hs*eno2(f0.tau,f1.tau,f2.tau,f3.tau),
		hs*eno2(f0.chi,f1.chi,f2.chi,f3.chi));
}

/** Calculates the ENO derivative using a sequence of values at
 * four gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double bd_sim::eno2(double p0,double p1,double p2,double p3) {
	return fabs(p0-2*p1+p2)>fabs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void bd_sim::write_files(int k) {
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
/*	if(fflags&2048) {
		double Q,max_qs;
		qs_measure(Q,max_qs);
		printf("QS info %g %.12g %.12g\n",time,Q,max_qs);
		output("qs",11,k);
	}
	if(fflags&12288) {
		compute_strain();
		if(fflags&4096) output("VS",11,k);
		if(fflags&8192) output("DS",12,k);
	}
	if(fflags&114688) {
		compute_def_rate();
		if(fflags&16384) output("Dtot",11,k);
		if(fflags&32768) output("Dpl",12,k);
		if(fflags&65536) output("Del",13,k);
	}*/
	if(fflags&131072) output("cc",14,k);
	if(fflags&262144) output("d",15,k);
	if(fflags&524288) output("c",16,k);
	if(fflags&1048576) output("phi",17,k);
	if(fflags&(1<<21)) print_extrema(k);
	if(eforce!=NULL) eforce->write(k);
}

/** Prints the extrema to a file. */
void bd_sim::print_extrema(int k) {

	// Open the output file. If this is the first frame, then remove
	// anything that was previously there. Otherwise append to the file.
	char *bufc=((char*) buf);
	sprintf(bufc,"%s/ex_file",filename);
	FILE *outf=safe_fopen(bufc,k==0?"w":"a");

	// Compute the extrema
	double gex[12];
	extrema(gex);

	double xmax=-100,ymax=-100;
	for(int j=0;j<n;j++) for(int i=0;i<m;i++) if(c[i+m*j]<4) {
		double x=ax+i*dx,y=ay+j*dy;
		if(x>xmax) xmax=x;
		if(y>ymax) ymax=y;
	}

	// Print the extrema and close the output file
	fprintf(outf,"%d %.14g %.14g %.14g %.14g %.14g %.14g %.14g "
		     "%.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g\n",
		     k,time,time*t_scale,*gex,gex[1],gex[2],gex[3],gex[4],
		     gex[5],gex[6],gex[7],gex[8],gex[9],gex[10],gex[11],xmax,ymax);
	fclose(outf);
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void bd_sim::output(const char *prefix,const int mode,const int sn) {

	// Determine whether to output a staggered field or not
	bool st=211452&(1<<mode);
	int *op=st?obnd+4:obnd;
	int i,j,l=op[1]-*op;

	// Assemble the output filename and open the output file
	char *bufc=((char*) buf);
	sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
	FILE *outf=safe_fopen(bufc,"wb");

	// Output the first line of the file
	float *bp=buf+1,*be=bp+l;
	*buf=l;
	for(i=*op;i<op[1];i++) *(bp++)=ax+(st?(i+0.5)*dx:i*dx);
	fwrite(buf,sizeof(float),l+1,outf);

	// Output the field values to the file
	if(mode<14) {

		// Deal with the regular fields using a common field pointer
		c_field *fp;
		for(j=op[2];j<op[3];j++) {
			fp=fm+me*j+*op;
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
				case 11: while(bp<be) *(bp++)=(fp++)->cu;break;
				case 12: while(bp<be) *(bp++)=(fp++)->cv;break;
				case 13: while(bp<be) *(bp++)=(fp++)->cp;
			}
			fwrite(buf,sizeof(float),l+1,outf);
		}
	} else {

		// Separately handle the fields relating to the level set
		for(j=op[2];j<op[3];j++) {
			*buf=ay+(st?(j+0.5)*dy:j*dy);bp=buf+1;
			switch(mode) {
				case 14: {int *ccp=cc+me*j+*op;while(bp<be) *(bp++)=*(ccp++);break;}
				case 15: {int *dp=qsm->d+me*j+*op;while(bp<be) *(bp++)=*(dp++);break;}
				case 16: {int *cp=c+m*j+*op;while(bp<be) *(bp++)=*(cp++);break;}
				case 17: {double *phip=phi+m*j+*op;while(bp<be) *(bp++)=*(phip++);}
			}
			fwrite(buf,sizeof(float),l+1,outf);
		}
	}

	// Close the file
	fclose(outf);
}

/** Sets the output routine to only output a portion of the simulation fields.
 * \param[in] (xmin,xmax) the x range to output.
 * \param[in] (ymin,ymax) the y range to output. */
void bd_sim::set_subfield_bounds(double xmin,double xmax,double ymin,double ymax) {
	double z;
	int q;

	// Set x subfield bounds and check that they are valid
	q=int(z=(xmin-ax)*xsp);*obnd=q<0?0:q;
	q=int(z-0.5);obnd[4]=q<0?0:q;
	q=int(z=(xmax-ax)*xsp+2);obnd[1]=q>me?me:q;
	q=int(z-0.5);obnd[5]=q>m?m:q;
	if(obnd[1]-*obnd<=0) fatal_error("Regular output grid has zero or negative width",1);
	if(obnd[5]-obnd[4]<=0) fatal_error("Staggered output grid has zero or negative width",1);

	// Set y subfield bounds and check that they are valid
	q=int(z=(ymin-ay)*ysp);obnd[2]=q<0?0:q;
	q=int(z-0.5);obnd[6]=q<0?0:q;
	q=int(z=(ymax-ay)*ysp+2);obnd[3]=q>ne?ne:q;
	q=int(z-0.5);obnd[7]=q>n?n:q;
	if(obnd[3]-obnd[2]<=0) fatal_error("Regular output grid has zero or negative height",1);
	if(obnd[7]-obnd[6]<=0) fatal_error("Staggered output grid has zero or negative height",1);
}

/** Computes the measure of the degree of quasi-staticity.
 * \param[out] Q the measure, integrated over the entire domain.
 * \param[out] max_qs the maximum value of the integrand.
 * \param[in] staggered whether to use the staggered velocity, only appropriate
 *			for quasi-static simulations. */
void bd_sim::qs_measure(double &Q,double &max_qs,bool staggered) {
	const double qsf1=6*K*mu*(K+(1/3.0)*mu),qsf2=6*K*mu*mu;
	Q=max_qs=0;

	// Deal with the case when the projection setup class hasn't been
	// allocated yet
	int qm;
	if(qsm==NULL) {
		for(c_field *fp=fm;fp<fm+mne;fp++) fp->cu=0;
		Q=max_qs=0;
		return;
	} else qm=qsm->m;

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		double ux,uy,vx,vy,l=0,max_l=0,l1,l2,l3;

		// Evaluate the gradient of velocity
		for(int i=0;i<m;i++) {
			c_field *fp=fm+i+me*j;

			// Skip any point outside of the current domain or
			// where the level set field is positive
			if(j<dom_lj||j>=dom_uj||i<dom_li||i>=dom_ui||c[i+m*j]>=4) {
				fp->cu=0;
				continue;
			}

			if(staggered) {
				vec *vp=vel+(j-dom_lj+m_pad)*qm+(i-dom_li);
				ux=0.5*xsp*(-vp->x+vp[1].x-vp[qm].x+vp[qm+1].x),
				uy=0.5*ysp*(-vp->x-vp[1].x+vp[qm].x+vp[qm+1].x);
				vx=0.5*xsp*(-vp->y+vp[1].y-vp[qm].y+vp[qm+1].y);
				vy=0.5*ysp*(-vp->y-vp[1].y+vp[qm].y+vp[qm+1].y);
			} else {
				c_field &fr=*(fp+1),&fl=*(fp-1);
				c_field &fd=*(fp-me),&fu=*(fp+me);
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

/** Computes the magnitudes of the deformation rate tensors and stores them in the  */
void bd_sim::compute_def_rate() {
	const double small_dev_cutoff=1e-12;
	double dt=dx*dx*tmult,itdtmu=1/(2*dt*mu),it_scale=1/t_scale;

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		int *cp=c+m*j;
		mat F,FT,C;
		for(c_field *fp=fm+me*j,*fe=fp+m;fp<fe;fp++,cp++) {
			c_field &f=*fp;
			if(*cp<4) {

				double ux=0.5*xsp*(fp[me+1].u-fp[me].u+fp[1].u-f.u),
				       uy=0.5*ysp*(fp[me+1].u+fp[me].u-fp[1].u-f.u),
				       vx=0.5*xsp*(fp[me+1].v-fp[me].v+fp[1].v-f.v),
				       vy=0.5*ysp*(fp[me+1].v+fp[me].v-fp[1].v-f.v),
				       t1,t2,t3,t4,t5,chit;

				// Total deformation rate
				f.cu=sqrt(0.5*(ux*ux+vy*vy)+0.25*(uy+vx)*(uy+vx))*it_scale;

				// Plastic deformation rate
				chit=f.chi;
				adaptive_plastic_term(f.dev(),chit,t4,dt);
				t4*=itdtmu;
				f.cv=t4*it_scale;

				// Elastic deformation rate
				t3=f.dev();
				t4/=(t3<small_dev_cutoff?small_dev_cutoff:t3);
				t1=ux-(f.s-f.q)*t4;
				t2=0.5*(uy+vx)-f.tau*t4;
				t3=vy-(-f.s-f.q)*t4;
				t5=2*f.q*t4;
				f.cp=sqrt(0.5*(t1*t1+t3*t3+t5*t5)+t2*t2)*it_scale;
			} else f.cu=f.cv=f.cp=0;
		}
	}
}

/* Computes invariants of the strain tensor. */
void bd_sim::compute_strain() {

	// Extrapolate the reference map
	ls.extrapolate_staggered_fields(e_rm);

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		int *cp=c+m*j;
		mat F,FT,C;
		for(c_field *fp=fm+me*j,*fe=fp+m;fp<fe;fp++,cp++) if(*cp<4) {
			double Xx=0.5*xsp*(fp[me+1].X-fp[me].X+fp[1].X-fp->X),
			       Xy=0.5*ysp*(fp[me+1].X+fp[me].X-fp[1].X-fp->X),
			       Yx=0.5*xsp*(fp[me+1].Y-fp[me].Y+fp[1].Y-fp->Y),
			       Yy=0.5*ysp*(fp[me+1].Y+fp[me].Y-fp[1].Y-fp->Y),
			       detF=1/(Xx*Yy-Xy*Yx),tr;

			// Compute the Green--Saint-Venant strain tensor
			F=mat(-Xx*detF,Yx*detF,Xy*detF,-Yy*detF);
			FT=F.transpose();
			C=(FT*F)*0.5;C.a-=0.5;C.d-=0.5;

			// Calculate and store invariants of the tensor
			tr=C.trace();
			C.a-=0.5*tr;C.d-=0.5*tr;
			fp->cu=detF;fp->cv=C.mod_sq();
		} else fp->cu=fp->cv=0;
	}
}

/** Calculates the extremal values of pressure, deviatoric stress, and effective
 * temperature.
 * \param[in] gex a pointer to 12-element array in which to store the values.
 *		  The first six are minima and maxima of the smoothed fields
 *		  (to remove any checkerboarding), and the next six are minima
 *		  and maxima of the original fields. */
void bd_sim::extrema(double *gex) {
	extrema_init(gex);

	// Find any gridpoint within the body
#pragma omp parallel for
	for(int j=0;j<m;j++) {
		double lex[12],*lp,de,ps,des,chis;
		c_field *fp=fm+me*j;
		int *cp=c+m*j;
		extrema_init(lex);
		for(int i=0;i<m;i++) if(cp[i]<4) {
			c_field *fpp=fp+i;

			// Compute smoothed values of the fields
			de=fpp->devsq();
			if(i>0&&cp[i-1]<4&&i<m-1&&cp[i+1]<4&&j>0&&cp[i-m]<4&&j<m-1&&cp[i+m]<4) {
				ps=fpp[-me].p+fpp[-1].p+fpp[1].p+fpp[me].p+4*fpp->p;
				des=fpp[-me].devsq()+fpp[-1].devsq()+fpp[1].devsq()+fpp[me].devsq()+4*de;
				chis=fpp[-me].chi+fpp[-1].chi+fpp[1].chi+fpp[me].chi+4*fpp->chi;
				lp=lex;
				extrema_check(0.125*ps,lp);
				extrema_check(0.125*des,lp);
				extrema_check(0.125*chis,lp);
			} else lp=lex+6;

			// Update local extrema counters
			extrema_check(fpp->p,lp);
			extrema_check(de,lp);
			extrema_check(fpp->chi,lp);
		}

#pragma omp critical
		{
			double *p=lex,*pe=p+12,*pg=gex;
			while(p<pe) {
				if(*p<*pg) *pg=*p;
				p++;pg++;
				if(*p>*pg) *pg=*p;
				p++;pg++;
			}
		}
	}

	// Take square root of deviatoric stresses and convert effective
	// temperatures to Kelvin
	gex[2]=sqrt(gex[2]);gex[3]=sqrt(gex[3]);
	gex[8]=sqrt(gex[8]);gex[9]=sqrt(gex[9]);
	gex[4]*=TZ;gex[5]*=TZ;
	gex[10]*=TZ;gex[11]*=TZ;
}
