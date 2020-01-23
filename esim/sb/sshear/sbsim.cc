#include "sbsim.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *		      vertical directions.
 * \param[in] (ax_,bx_) the lower and upper x-coordinate simulation bounds.
 * \param[in] (ay_,by_) the lower and upper y-coordinate simulation bounds.
 * \param[in] chi_inf_ the parameter \f$\chi_\infty\f$ in the STZ model.
 * \param[in] c_0_inv_ the parameter \f$1/c_0\f$ in the STZ model.
 * \param[in] nu_ the parameter \f$\nu\f$ in the STZ model.
 * \param[in] mu_ the elastic shear modulus.
 * \param[in] K_ the elastic bulk modulus.
 * \param[in] visc_ the viscosity.
 * \param[in] tmult_ a multiplier to apply to the default timestep computation.
 * \param[in] filename_ the filename of the output directory. */
sbsim::sbsim(int m_,int n_,double ax_,double bx_,double ay_,double by_,
	     double chi_inf_,double c_0_inv_,double nu_,double mu_,
	     double K_,const double visc_,double tmult_,const char *filename_)
	: m(m_), n(n_), mn(m_*n_), ax(ax_), ay(ay_), dx((bx_-ax_)/m_),
	dy((by_-ay_)/n_), xsp(1/dx), ysp(1/dy), chi_inf(chi_inf_),
	c_0_inv(c_0_inv_), nu(nu_), mu(mu_), K(K_), mu_inv(1/mu),
	tmult(tmult_), visc(visc_), filename(filename_) {

	// Dynamically allocate memory for the simulation fields
	u=new double[mn];v=new double[mn];p=new double[mn];
	s=new double[mn],tau=new double[mn];chi=new double[mn];

	// Dynamically allocate memory for the field updates
	cu=new double[mn];cv=new double[mn];cp=new double[mn];
	cs=new double[mn];ctau=new double[mn];cchi=new double[mn];

	// Initialize the array used for output
	buf=new float[m+1];bufe=buf+m+1;

	// Initialize the simulation fields
	init_fields();
}

/** The class destructor frees the dynamically allocated memory. */
sbsim::~sbsim() {
	delete [] buf;
	delete [] cchi;delete [] ctau;delete [] cs;
	delete [] cp;delete [] cv;delete [] cu;
	delete [] chi;delete [] tau;delete [] s;
	delete [] p;delete [] v;delete [] u;
}

/** Carries out the simulation for a specified time interval,
 * periodically saving the output.
 * \param[in] t_start the simulation time to start from.
 * \param[in] t_end the simulation time to end at.
 * \param[in] frames the number of frames to save. */
void sbsim::solve(double t_start,double t_end,int frames) {
	time=t_start;
	double time_interval=(t_end-t_start)/frames,target_time,dt=dx*dx*tmult;

	// Output the initial fields
	write_files(0);
	puts("# Output frame 0");

	for(int k=1;k<=frames;k++) {

		// Compute the target time to the next output frame
		target_time=t_start+time_interval*k;

		// Carry out simulation step using the regular timestep until
		// within range of the target time
		while(time+dt+1e-11<target_time) step_forward(dt);

		// Carry out a final simulation step, using exactly the right
		// time step to reach the target time
		step_forward(target_time-time);

		// Output the fields
		write_files(k);
		printf("# Output frame %d\n",k);
	}

}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void sbsim::step_forward(double dt) {
	double hx,hy,x,y;
	double u1,v1,p1,s1,tau1,chi1;
	double u2,v2,p2,s2,tau2,chi2;
	double ux,vx,px,sx,taux;
	double uy,vy,py,sy,tauy;
	double tomega,magn,magr,cha;
	double hxx,hyy,uxx,vxx;
	int i,j,ij,ijd,iju,ijdd,ijuu;

	// Update the simulation time
	time+=dt;

	// Compute grid spacing constants based on the current timestep
	hx=0.5*dt*xsp;hy=0.5*dt*ysp;
	hxx=2*hx*xsp;hyy=2*hy*ysp;

	// Loop over the grid points
	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) {

		// Find the locations of the grid points to the left and right
		// of this one, taking into account periodicity
		ijd=i==0?ij+m-1:ij-1;
		iju=i==m-1?ij-m+1:ij+1;

		// Calculate a one-sided derivative using ENO2 if the
		// gridpoints are available
		if(u[ij]>0) {
			ijdd=i<=1?ij+m-2:ij-2;
			eno2(iju,ij,ijd,ijdd,u1,v1,p1,s1,tau1,chi1);
			u1*=hx;v1*=hx;p1*=hx;s1*=hx;tau1*=hx;chi1*=hx;
		} else {
			ijuu=i>=m-1?ij-m+2:ij+2;
			eno2(ijd,ij,iju,ijuu,u1,v1,p1,s1,tau1,chi1);
			u1*=-hx;v1*=-hx;p1*=-hx;s1*=-hx;tau1*=-hx;chi1*=-hx;
		}

		// Compute first-order centered differences
		tomega=hx*(v[iju]-v[ijd]);
		ux=hx*(u[iju]-u[ijd]);
		vx=hx*(v[iju]-v[ijd]);
		px=hx*(p[iju]-p[ijd]);
		sx=hx*(s[iju]-s[ijd]);
		taux=hx*(tau[iju]-tau[ijd]);

		// Compute viscosity contribution
		uxx=hxx*(u[iju]-2*u[ij]+u[ijd]);
		vxx=hxx*(v[iju]-2*v[ij]+v[ijd]);

		// Find the locations of the grid points that are up and down
		// of this one, taking into account periodicity
		ijd=j==0?ij+mn-m:ij-m;
		iju=j==n-1?ij-mn+m:ij+m;

		// Calculate a one-sided derivative using ENO2
		// if the gridpoints are available
		if(v[ij]>0) {
			ijdd=j<=1?ij+mn-2*m:ij-2*m;
			eno2(iju,ij,ijd,ijdd,u2,v2,p2,s2,tau2,chi2);
			u2*=hy;v2*=hy;p2*=hy;s2*=hy;tau2*=hy;chi2*=hy;
		} else {
			ijuu=j>=n-2?ij-mn+2*m:ij+2*m;
			eno2(ijd,ij,iju,ijuu,u2,v2,p2,s2,tau2,chi2);
			u2*=-hy;v2*=-hy;p2*=-hy;s2*=-hy;tau2*=-hy;chi2*=-hy;
		}

		// Compute first order centered-differences
		tomega+=hy*(u[ijd]-u[iju]);
		uy=hy*(u[iju]-u[ijd]);
		vy=hy*(v[iju]-v[ijd]);
		py=hy*(p[iju]-p[ijd]);
		sy=hy*(s[iju]-s[ijd]);
		tauy=hy*(tau[iju]-tau[ijd]);

		// Compute update to velocity and pressure
		cu[ij]=-u[ij]*u1-v[ij]*u2+(-px+sx+tauy)*mu_inv;
		cv[ij]=-u[ij]*v1-v[ij]*v2+(-py-sy+taux)*mu_inv;
		cp[ij]=-u[ij]*p1-v[ij]*p2-K*(ux+vy);

		// Compute the STZ plasticity and the update to chi
		magn=s[ij]*s[ij]+tau[ij]*tau[ij];
		magr=sqrt(magn);
		cha=dt*nu*exp(-1/chi[ij])*qs(magr);
		cchi[ij]=(chi[ij]-u[ij]*chi1-v[ij]*chi2+cha*chi_inf*c_0_inv)/(1+cha*c_0_inv);
		magr=1/(1+mu*cha/(std::abs(magr)>1e-8?magn:1));

		// Compute the update to the deviatoric stress fields
		cs[ij]=(s[ij]-u[ij]*s1-v[ij]*s2-tomega*tau[ij]+mu*(ux-vy))*magr;
		ctau[ij]=(tau[ij]-u[ij]*tau1-v[ij]*tau2+tomega*s[ij]+mu*(uy+vx))*magr;

		// Add the viscous contribution
		cu[ij]+=visc*(uxx+hyy*(u[iju]-2*u[ij]+u[ijd]));
		cv[ij]+=visc*(vxx+hyy*(v[iju]-2*v[ij]+v[ijd]));

		// Add in the Gaussian forcing
		x=ax+(i+0.5)*dx;
		y=ay+(j+0.5)*dy;
		cv[ij]-=dt*0.5*exp(-8*((x-0.5)*(x-0.5)+y*y));
		cv[ij]+=dt*0.5*exp(-8*((x+0.5)*(x+0.5)+y*y));
	}

	// Apply the changes
	for(ij=0;ij<mn;ij++) {
		u[ij]+=cu[ij];cu[ij]=0;
		v[ij]+=cv[ij];cv[ij]=0;
		p[ij]+=cp[ij];cp[ij]=0;
		s[ij]=cs[ij];cs[ij]=0;
		tau[ij]=ctau[ij];ctau[ij]=0;
		chi[ij]=cchi[ij];cchi[ij]=0;
	}
}

/** Computes and prints the extrema of several simulation fields. */
void sbsim::print_extrema() {
	double minu,maxu,minp,maxp,minchi,maxchi;
	minu=maxu=u[0];
	minp=maxp=p[0];
	minchi=maxchi=chi[0];
	for(int ij=1;ij<mn;ij++) {
		if(u[ij]<minu) minu=u[ij];
		if(u[ij]>maxu) maxu=u[ij];
		if(p[ij]<minp) minp=p[ij];
		if(p[ij]>maxp) maxp=p[ij];
		if(chi[ij]<minchi) minchi=chi[ij];
		if(chi[ij]>maxchi) maxchi=chi[ij];
	}
	printf("Ext: %f %f %f %f %f %f %f\n",time,minu,maxu,minp,maxp,minchi,maxchi);
}

/** Computes the ENO2 derivative for all the simulation fields.
 * \param[in] (ij0,ij1,ij2,ij3) the four grid points to consider.
 * \param[out] (uv,vv,pv,sv,tauv,chiv) the computed derivatives, without any
 *                                     normalization for the grid spacing. */
void sbsim::eno2(int ij0,int ij1,int ij2,int ij3,double &uv,double &vv,double &pv,double &sv,double &tauv,double &chiv) {
	uv=eno2(u[ij0],u[ij1],u[ij2],u[ij3]);
	vv=eno2(v[ij0],v[ij1],v[ij2],v[ij3]);
	pv=eno2(p[ij0],p[ij1],p[ij2],p[ij3]);
	sv=eno2(s[ij0],s[ij1],s[ij2],s[ij3]);
	tauv=eno2(tau[ij0],tau[ij1],tau[ij2],tau[ij3]);
	chiv=eno2(chi[ij0],chi[ij1],chi[ij2],chi[ij3]);
}

/** Computes the ENO2 finite difference scheme for four values.
 * \param[in] (p0,p1,p2,p3) the four values to use.
 * \return The computed ENO2 derivative. */
inline double sbsim::eno2(double p0,double p1,double p2,double p3) {
	return std::abs(p0-2*p1+p2)>std::abs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

/** Computes the quantity s*q(s) in the athermal STZ model.
 * \param[in] s the input magnitude of shear stress.
 * */
inline double sbsim::qs(double s) {
	if(s<0) s=-s;
	return s>1?(s-1)*(s-1):0;
}

/** Writes a selection of simulation fields to the output directory.
 * \param[in] k the frame number to append to the output. */
void sbsim::write_files(int k) {
	output("u",u,k);
	output("v",v,k);
	output("p",p,k);
	output("s",s,k);
	output("tau",tau,k);
	output("chi",chi,k);
	for(int ij=0;ij<mn;ij++) cu[ij]=sqrt(s[ij]*s[ij]+tau[ij]*tau[ij]);
	output("dev",cu,k);
}

/** Saves a simulation field to the output directory in a binary format that
 * can be read by Gnuplot.
 * \param[in] fname the filename to use.
 * \param[in] fp a pointer to the field to save.
 * \param[in] sn the frame number to append to the filename. */
void sbsim::output(const char *fname,double *fp,const int sn) {
	int i,j;
	float *bp;

	// Create the output filename and open the file
	char *fn=new char[strlen(fname)+20];
	sprintf(fn,"%s/%s.%d",filename,fname,sn);
	FILE *outf=safe_fopen(fn,"wb");
	delete [] fn;

	// Output the header line
	*buf=m;bp=buf+1;
	for(i=0;i<m;i++) *(bp++)=ax+(i+0.5)*dx;
	fwrite(buf,sizeof(float),m+1,outf);

	// Output each line of the field, converting from double precision to
	// single precision
	for(j=0;j<n;j++) {
		*buf=ay+(j+0.5)*dy;bp=buf+1;
		while(bp<bufe) *(bp++)=*(fp++);
		fwrite(buf,sizeof(float),m+1,outf);
	}

	fclose(outf);
}
