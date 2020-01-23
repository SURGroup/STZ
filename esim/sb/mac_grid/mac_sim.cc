#include <cstring>

#include "common.hh"
#include "mac_sim.hh"

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
 * \param[in] tmult_ a multiplier to apply to the default timestep size.
 * \param[in] filename_ the filename of the output directory. */
mac_sim::mac_sim(const int m_,const int n_,const double ax_,const double bx_,
		const double ay_,const double by_,const double mu_,
		const double K_,const double viscosity_,
		const double tmult_,const char *filename_)
	: m(m_), n(n_), mn(m_*n_), ml(m+4), ax(ax_), ay(ay_), bx(bx_), by(by_),
	dx((bx_-ax_)/m_), dy((by_-ay_)/n_), xsp(1/dx), ysp(1/dy), mu(mu_),
	mu_inv(1/mu_), K(K_), viscosity(viscosity_), tmult(tmult_),
	filename(filename_), fbase(new field[ml*(n+4)]), fm(fbase+2*ml+2),
	time(0), f_num(0), buf(new float[m>123?m+5:128]) {}

/** The class destructor frees the dynamically allocated memory. */
mac_sim::~mac_sim() {
	delete [] buf;
	delete [] fbase;
}

/** Initializes the simulation fields. */
void mac_sim::init_fields() {

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		double y=ay+dy*(j+0.5),x=ax+dx*0.5;
		for(field *fp=fm+ml*j,*fe=fp+m;fp<fe;fp++) {
			fp->u=1.5*exp(-16*(x*x+y*y));
			fp->v=0;
			fp->clear_stress();
			fp->X=x;
			fp->Y=y;
			x+=dx;
		}
	}

	// Now that the primary grid points are set up, initialize the ghost
	// points according to the boundary conditions
	set_boundaries();
}

/** Carries out the simulation for a specified time interval using the direct
 * simulation method, periodically saving the output.
 * \param[in] duration the simulation time to end at.
 * \param[in] frames the number of frames to save. */
void mac_sim::solve(double duration,int frames) {
	const double dt=dx*dx*tmult;
	double t_start=time,time_interval=duration/frames;
	int l=int(time_interval/dt)+1;
	double adt=time_interval/l,t0,t1,t2;

	// Output the initial fields and record initial time
	if(f_num==0) {write_files(0);puts("# Output frame 0");}
	t0=wtime();

	for(int k=1;k<=frames;k++) {

		// Perform the required number of timesteps
		for(int i=0;i<l;i++) step_forward(adt);

		// Output the fields
		t1=wtime();
		write_files(k+f_num);

		// Print diagnostic information, and set the time manually to
		// correct for any minor flaoting point errors
		time=t_start+k*time_interval;
		t2=wtime();
		printf("# Output frame %d [%d, %.8g s, %.8g s]\n",k+f_num,l,t1-t0,t2-t1);
		t0=t2;
	}
	f_num+=frames;
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void mac_sim::step_forward(double dt) {
	const double lambda=(K-(2/3.)*mu),rho_inv=1;
	double hx=0.5*xsp*dt,hy=0.5*ysp*dt;
	double hxxv=2*xsp*hx*viscosity,hyyv=2*ysp*hy*viscosity,hcv=2*(hxxv+hyyv);

	// Update the simulation time
	time+=dt;

#pragma omp parallel for
	for(int j=0;j<n;j++) {
		for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++) {

			// Create references to the fields in the neighboring gridpoints,
			// taking into account the periodicity
			field &f=*fp,&fl=fp[-1],&fr=fp[1],&fd=fp[-ml],&fu=fp[ml];
			double uc=f.u,vc=f.v,ul=0.5*(f.u+fl.u),vl=0.5*(f.v+fl.v),
			       ud=0.5*(f.u+fd.u),vd=0.5*(f.v+fd.v),
			       ux,vx,uy,vy,s11x,s11y,s12x,s12y,s21x,s21y,s22x,s22y,
			       s33x,s33y,Xx,Xy,Yx,Yy,uvisc,vvisc,fx,fy,dxy;

			// Calculate ENO derivatives of the fields
			uc>0?c_eno2(ux,vx,s33x,Xx,Yx,hx,fr,f,fl,fp[-2])
			    :c_eno2(ux,vx,s33x,Xx,Yx,-hx,fl,f,fr,fp[2]);
			vc>0?c_eno2(uy,vy,s33y,Xy,Yy,hy,fu,f,fd,fp[-2*ml])
			    :c_eno2(uy,vy,s33y,Xy,Yy,-hy,fd,f,fu,fp[2*ml]);

			// Calculate net force due to stress imbalance
			net_force(fp,hx,hy,fx,fy);

			// Calculate viscosity
			uvisc=hxxv*(fr.u+fl.u)+hyyv*(fu.u+fd.u)-hcv*uc;
			vvisc=hxxv*(fr.v+fl.v)+hyyv*(fu.v+fd.v)-hcv*vc;

			// Calculate the update to velocity
			f.cu=-uc*ux-vc*uy+rho_inv*(fx+uvisc);
			f.cv=-uc*vx-vc*vy+rho_inv*(fy+vvisc);

		 	// Calculate update to reference map fields
			f.cX=-uc*Xx-vc*Xy;
			f.cY=-uc*Yx-vc*Yy;

			// Calculate updates to out-of-plane stress component
			ux=hx*(fp[1].u-fp[-1].u);
			vy=hy*(fp[ml].v-fp[-ml].v);
			f.cs33=-uc*s33x-vc*s33y+lambda*(ux+vy);

			// Compute advective terms on the left edge
			ul>0?l_eno2(s11x,s21x,hx,fr,f,fl,fp[-2])
			    :l_eno2(s11x,s21x,-hx,fl,f,fr,fp[2]);
			vl>0?l_eno2(s11y,s21y,hy,fu,f,fd,fp[-2*ml])
			    :l_eno2(s11y,s21y,-hy,fd,f,fu,fp[2*ml]);

			// Compute centered differences on the left edge
			ux=hx*(fp->u-fp[-1].u);
			vx=hx*(fp->v-fp[-1].v);
			uy=0.25*hy*(fp[ml].u+fp[ml-1].u-fp[-ml].u-fp[-ml-1].u);
			vy=0.25*hy*(fp[ml].v+fp[ml-1].v-fp[-ml].v-fp[-ml-1].v);

			// Update stresses on the left edge
			dxy=0.5*(vx+uy);
			f.cs11=-ul*s11x-vl*s11y+lambda*(ux+vy)+2*mu*ux;
			f.cs21=-ul*s21x-vl*s21y+2*mu*dxy;

			// Compute advective terms on the bottom edge
			ud>0?d_eno2(s12x,s22x,hx,fr,f,fl,fp[-2])
			    :d_eno2(s12x,s22x,-hx,fl,f,fr,fp[2]);
			vd>0?d_eno2(s12y,s22y,hy,fu,f,fd,fp[-2*ml])
			    :d_eno2(s12y,s22y,-hy,fd,f,fu,fp[2*ml]);

			// Compute centered differences on the bottom edge
			ux=0.25*hx*(fp[1-ml].u+fp[1].u-fp[-1].u-fp[-1-ml].u);
			vx=0.25*hx*(fp[1-ml].v+fp[1].v-fp[-1].v-fp[-1-ml].v);
			uy=hy*(fp->u-fp[-ml].u);
			vy=hy*(fp->v-fp[-ml].v);

			// Update stresses on the bottom edge
			dxy=0.5*(vx+uy);
			f.cs12=-ud*s12x-vd*s12y+2*mu*dxy;
			f.cs22=-ud*s22x-vd*s22y+lambda*(ux+vy)+2*mu*vy;
		}
	}

	// Apply the updates to the fields
#pragma omp parallel for
	for(int j=0;j<n;j++) {
		for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++) fp->update();
	}

	// Set the top and bottom boundaries
	set_boundaries();
}

/** Calculates the net force at a grid point. */
inline void mac_sim::net_force(field *fp,double xf,double yf,double &fx,double &fy) {
	field &f=*fp,&fr=fp[1],&fu=fp[ml];
	fx=(fr.s11-f.s11)*xf+(fu.s12-f.s12)*yf;
	fy=(fr.s21-f.s21)*xf+(fu.s22-f.s22)*yf;
}

/** Sets the fields in the ghost regions according to the boundary conditions.
 */
void mac_sim::set_boundaries() {
	const double lx=bx-ax,ly=by-ay;

	// Set left and right ghost values
	for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
		fp[-2].copy(fp[m-2],-lx,0);
		fp[-1].copy(fp[m-1],-lx,0);
		fp[m].copy(*fp,lx,0);
		fp[m+1].copy(fp[1],lx,0);
	}

	// Set top and bottom ghost values
	const int tl=2*ml,g=m*ml;
	for(field *fp=fm-2,*fe=fp+ml;fp<fe;fp++) {
		fp[-tl].copy(fp[g-tl],0,-ly);
		fp[-ml].copy(fp[g-ml],0,-ly);
		fp[g].copy(*fp,0,ly);
		fp[g+ml].copy(fp[ml],0,ly);
	}
}

/** Calculates one-sided derivatives of the fields using the second-order ENO2
 * scheme, applying the shift to the X terms.
 * \param[in] hs a multiplier to apply to the computed fields.
 * \param[in] (f0,f1,f2,f3) the fields to compute the derivative with.
 * \param[in] (X0,X2,X3) corresponding shifts to apply to the X terms.
 * \return The computed derivative. */
inline void mac_sim::c_eno2(double &uv,double &vv,double &s33v,double &Xv,double &Yv,double hs,field &f0,field &f1,field &f2,field &f3) {
	uv=hs*eno2(f0.u,f1.u,f2.u,f3.u);
	vv=hs*eno2(f0.v,f1.v,f2.v,f3.v);
	s33v=hs*eno2(f0.s33,f1.s33,f2.s33,f3.s33);
	Xv=hs*eno2(f0.X,f1.X,f2.X,f3.X);
	Yv=hs*eno2(f0.Y,f1.Y,f2.Y,f3.Y);
}

inline void mac_sim::l_eno2(double &s11v,double &s21v,double hs,field &f0,field &f1,field &f2,field &f3) {
	s11v=hs*eno2(f0.s11,f1.s11,f2.s11,f3.s11);
	s21v=hs*eno2(f0.s21,f1.s21,f2.s21,f3.s21);
}

inline void mac_sim::d_eno2(double &s12v,double &s22v,double hs,field &f0,field &f1,field &f2,field &f3) {
	s12v=hs*eno2(f0.s12,f1.s12,f2.s12,f3.s12);
	s22v=hs*eno2(f0.s22,f1.s22,f2.s22,f3.s22);
}

/** Calculates the ENO derivative using a sequence of values at
 * four gridpoints.
 * \param[in] (p0,p1,p2,p3) the sequence of values to use.
 * \return The computed derivative. */
inline double mac_sim::eno2(double p0,double p1,double p2,double p3) {
	return abs(p0-2*p1+p2)>abs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void mac_sim::write_files(int k) {
	const int fflags=63;
	if(fflags&1) output("u",0,k);
	if(fflags&2) output("v",1,k);
	if(fflags&4) output("p",2,k,false);
	if(fflags&8) output("dev",3,k,false);
	if(fflags&16) output("X",4,k);
	if(fflags&32) output("Y",5,k);
}

/** Outputs a 2D array to a file in a format that can be read by Gnuplot.
 * \param[in] prefix the field name to use as the filename prefix.
 * \param[in] mode the code of the field to print.
 * \param[in] sn the current frame number to append to the filename. */
void mac_sim::output(const char *prefix,const int mode,const int sn,const bool ghost) {

	// Determine whether to output a cell-centered field or not
	bool cen=true;
	int l=ghost?ml:(cen?m:m+1);
	double disp=(cen?0.5:0)-(ghost?2:0);

	// Assemble the output filename and open the output file
	char *bufc=((char*) buf);
	sprintf(bufc,"%s/%s.%d",filename,prefix,sn);
	FILE *outf=safe_fopen(bufc,"wb");

	// Output the first line of the file
	int i,j;
	float *bp=buf+1,*be=bp+l;
	*buf=l;
	for(i=0;i<l;i++) *(bp++)=ax+(i+disp)*dx;
	fwrite(buf,sizeof(float),l+1,outf);

	// Output the field values to the file
	field *fr=ghost?fbase:fm;
	for(j=0;j<l;j++,fr+=ml) {
		field *fp=fr;
		*buf=ay+(j+disp)*dy;bp=buf+1;
		switch(mode) {
			case 0: while(bp<be) *(bp++)=(fp++)->u;break;
			case 1: while(bp<be) *(bp++)=(fp++)->v;break;
			case 2: while(bp<be) *(bp++)=pressure(fp++);break;
			case 3: while(bp<be) *(bp++)=dev(fp++);break;
			case 4: while(bp<be) *(bp++)=(fp++)->X;break;
			case 5: while(bp<be) *(bp++)=(fp++)->Y;break;
			case 6: while(bp<be) *(bp++)=(fp++)->s11;break;
			case 7: while(bp<be) *(bp++)=(fp++)->s12;break;
			case 8: while(bp<be) *(bp++)=(fp++)->s21;break;
			case 9: while(bp<be) *(bp++)=(fp++)->s22;break;
			case 10: while(bp<be) *(bp++)=(fp++)->s33;
		}
		fwrite(buf,sizeof(float),l+1,outf);
	}

	// Close the file
	fclose(outf);
}

double mac_sim::pressure(field *fp) {
	return -(1/6.)*(2*fp->s33+fp->s11+fp->s22+fp[1].s11+fp[ml].s22);
}

double mac_sim::dev_sq(field *fp) {
	double s11=0.5*(fp->s11+fp[1].s11),s21=0.5*(fp->s21+fp[1].s21),
	       s12=0.5*(fp->s12+fp[ml].s12),s22=0.5*(fp->s22+fp[ml].s22),
	       s33=fp->s33,p=-(1/3.)*(s11+s22+s33);
	s11+=p;s22+=p;s33+=p;
	return 0.5*(s11*s11+s12*s12+s21*s21+s22*s22+s33*s33);
}
