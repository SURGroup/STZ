#include <sys/types.h>
#include <sys/stat.h>

#include "sbfrac.hh"

void sbfrac::solve() {
	double t_end=max_ki/(gamma*typical_ki),time_interval=t_end/frames,target_time,t0,t1,t2;
	int k;

	fs=new bd_sim(m,m,-sim_size,sim_size,-sim_size,sim_size,
			mu,K,viscosity,chi_len,filter_stress,t_scale,0.5,
			adapt_fac,stz,filename,fflags);
	fs->bc_smooth=bc_smooth;
	if(limit_output) fs->set_subfield_bounds(-output_size,output_size,-output_size,output_size);
	double dx=fs->dx,dt=tmult*dx*dx;
	init_fields();

	// Create output directory
	mkdir(filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Open diagnostic files
	sprintf(buf,"%s/kfile",filename);kfile=safe_fopen(buf,"wb");
	sprintf(buf,"%s.fld",filename);ffile=safe_fopen(buf,"wb");

	// Write initial frame
	fs->write_files(0);puts("# Output frame 0");
	t0=wtime();

	for(k=1;k<=frames;k++) {

		// Determine the integer number of frames required in order to
		// make sure the timestep is less than dt
		int l=int(time_interval/dt)+1;
		double adt=time_interval/l;

		// Perform the timesteps, checking for void nucleation
		for(int ll=0;ll<l;ll++) {
			step_forward_quasistatic(adt);
			if(nucleated) break;
		}
		if(nucleated) break;

		// Adjust timestep if adaptive option is selected
		if(adaptive_ts) {
			double dtn=adaptive_ts_select();
			if(dtn>0&&dt>dtn) dt=dtn;
		}
		t1=wtime();

		// Output files
		fs->write_files(k);
		t2=wtime();
		printf("# Output frame %d, t=%g [%d steps, com=%.8g s, write=%.8g s]\n",k,fs->time,l,t1-t0,t2-t1);
		t0=t2;
	}

	// If void nucleation has occurred, then carry out additional frames
	// using direct code
	if(nucleated) {
		dt=0.5*dx*dx;
		fs->ls.failsafe=true;
		printf("# Void nucleation output frame %d\n",k);
		fs->write_files(k);
		int kbase=k;double tbase=fs->time;
		for(k++;k<=kbase+frames;k++) {
			target_time=tbase+time_interval*void_nucl_tfac*(k-kbase);
			while(fs->time+dt*1.000001<target_time) step_forward(dt);
			step_forward(target_time-fs->time);

			// Output files
			printf("# Output frame %d\n",k);
			fs->write_files(k);
		}
	}

	// Close any open diagnostic files
	fclose(kfile);fclose(ffile);

	// Delete simulation class
	delete fs;
}

/** Calculates an adaptive timestep based on the maximum value of plastic
 * deformation. */
double sbfrac::adaptive_ts_select() {
	const double &dx=fs->dx;
	double dpex[16],&chi_len=fs->chi_len;
	fs->extrema(dpex);
	return dpex[15]>1e-30?dx*dx/(2*chi_len*chi_len*dpex[15]):bd_double_max;
}

/** Steps the simulation fields forward using the direct update procedure.
 * \param[in] dt the time step to use. */
void sbfrac::step_forward(double dt) {
	fs->time+=dt;
	fs->direct_step(dt);
	if(void_nucl) nucleate_void();
	post_process();
	if(++d_counter==100) {print_extrema();d_counter=0;}
}

/** Steps the simulation fields forward using the quasistatic update procedure.
 * \param[in] dt the time step to use. */
void sbfrac::step_forward_quasistatic(double dt) {

	// Update the simulation time
	fs->time+=dt;

	// Carry out the intermediate step that considers the advective and
	// plasticity terms
	for(int i=0;i<a_substep;i++) {
		fs->advection_step(dt/a_substep);
		post_process(false);
	}

	// Apply the projection step to ensure quasistaticity
	fs->projection_step(dt);
	if(void_nucl) nucleate_void();
	post_process();
	print_extrema();
}

/** Applies operations to the simulation fields after they have been updated by
 * time integration. It extrapolates the fields to neighboring gridpoints, and
 * it imposes boundary conditions at the edge of the simulation domain.
 * \param[in] all_fields whether to extrapolate all simulation fields, or just
 *			 the cell-centered ones. */
void sbfrac::post_process(bool all_fields) {
	fs->extrapolate(all_fields);
	set_boundaries();
}

/** Prints the extrema of the pressure, deviatoric stress and effective
 * temperature. */
void sbfrac::print_extrema() {
	double gex[16],dx=fs->dx,&chi_len=fs->chi_len;
	fs->extrema(gex);

	// Calculate the timesteps based on the maximum allowed from the
	// plastic deformation field. There are two variations: the first uses
	// the smoothed field, and the second uses the original field values.
	double dtn=(dx*dx)/(2*chi_len*chi_len*gex[7]),
	       dtn1=(dx*dx)/(2*chi_len*chi_len*gex[15]);

	// Output the results to the extrema file
	fprintf(ffile,"%.14g %.14g %.14g %.14g %.14g %.14g %.14g"
		      " %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g\n",
			fs->time,fs->time*t_scale,typical_ki*1e-6*(lambda+gamma*fs->time),
			*gex,gex[1],gex[2],gex[3],gex[4],gex[5],
			gex[6],gex[7],gex[8],gex[9],gex[10],gex[11],
			gex[12],gex[13],gex[14],gex[15],dtn,dtn1);
}

/** Searches for gridpoints with pressure values below a cutoff, and replaces
 * them with small voids. */
void sbfrac::nucleate_void() {
	bool any_voids=false;
	double *phi=fs->phi;
	c_field *fm=fs->fm;
	int *c=fs->c;

#pragma omp parallel for
	for(int j=0;j<m;j++) {
		for(int ij=j*m;ij<(j+1)*m;ij++) if(c[ij]<4&&fm[ij+j].p<void_nucl_p) {
#pragma omp critical
			{
				// Introduce the void
				int i=ij-j*m,ai=i>10?i-10:0,bi=i+11<m?i+11:m,
				    aj=j>10?j-10:0,bj=j+11<m?j+11:m,i2,j2;
				printf("# Void at (%d,%d)\n",i,j);
				double vox=i+rshift(),voy=j+rshift(),xs,ys,phit;
				for(j2=aj;j2<bj;j2++) for(i2=ai;i2<bi;i2++) {
					double &phic=phi[i2+j2*m];
					xs=vox-i2;ys=voy-j2;
					phit=fs->dx*(void_nucl_rad-sqrt(xs*xs+ys*ys));
					if(phit>phic) phic=phit;
					any_voids=true;
				}
			}
		}
	}

	// If any voids have been introduced, the level set narrow band must be
	// rebuilt
	if(any_voids) {nucleated=true;fs->ls.build_band();}
}

/** Writes a selection of output files.
 * \param[in] k the numerical index to append to the file output. */
void sbfrac::write_files(int k) {
	fs->write_files(k);

	// Save information about the frame number and stress intensity factor,
	// flushing the buffer to keep the file up-to-date
	fprintf(kfile,"%d %.14g %.14g %.14g\n",k,fs->time,fs->time*t_scale,typical_ki*1e-6*(lambda+gamma*fs->time));
	fflush(kfile);
}

/** Computes the argument of a given vector, breaking it up into several cases
 * to avoid doing arctangents of very large numbers.
 * \param[in] (x,y) the vector to consider.
 * \return the argument in radians, measured counterclockwise from the positive
 * x axis. */
double sbfrac::arg(double x,double y) {
	const double pi=3.1415926535897932384626433832795;
	return x+y>0?(x>y?atan(y/x):0.5*pi-atan(x/y)):(x>y?-atan(x/y)-0.5*pi:atan(y/x)+(y>0?pi:-pi));
}

void sbfrac::init_fields() {
	c_field *fm=fs->fm;
	int *c=fs->c,*cc=fs->cc,j,k;
	double fac,scale,*rnoise,*phip=fs->phi,dx=fs->dx;

	// Initialize the random chi field if needed
	if(random_chi0) random_chi_init(rnoise,k,fac,scale);

	// Set up the level set to describe the object boundary
	for(j=0;j<=m;j++) {
		double y=-sim_size+(j+0.5)*dx;
		for(int i=0;i<m;i++,phip++) {
			double x=-sim_size+(i+0.5)*dx;
			*phip=1-(x+displacement<0?fabs(y):sqrt((x+displacement)*(x+displacement)+y*y));

			// Create initial void if requested
			if(init_void) {
				x-=init_void_pos;
				double tphi=init_void_rad-sqrt(x*x+y*y);
				if(tphi>*phip) *phip=tphi;
			}

		}
	}

	// Build the band, which will set the status array, and use it to set
	// the staggered status array
	fs->ls.build_band();
	fs->set_staggered_status();

#pragma omp parallel for
	for(j=0;j<=m;j++) {
		double ys=-sim_size+j*dx,y=ys+0.5*dx;
		for(int i=0;i<=m;i++) {
			double xs=-sim_size+i*dx,x=xs+0.5*dx,r,theta;
			int ij=i+m*j,ije=i+me*j;
			c_field &f=fm[ije];

			// Set regular fields
			if(i<m&&j<m&&c[ij]<4) {
				r=sqrt(x*x+y*y);
				theta=r<1e-20?0:arg(x,y);

				// In-plane stress components
				r=Kinit/sqrt(r);
				f.p=mix_i*(-r*cos(theta*0.5))+mix_ii*(r*sin(theta*0.5));
				f.s=mix_i*(-r*cos(theta*0.5)*sin(theta*0.5)*sin(theta*1.5))
				   +mix_ii*(r*cos(theta*0.5)*sin(theta*0.5)*cos(theta*1.5));
				f.tau=mix_i*(r*sin(theta*0.5)*cos(theta*0.5)*cos(theta*1.5))
				     +mix_ii*(r*cos(theta*0.5)*(1-sin(theta*0.5)*sin(theta*1.5)));

				// Out-of-plane stress component and effective
				// temperature
				f.q=-poisson*f.p*0.25;
				f.chi=chi0;

				// Random chi modification if requested
				if(random_chi0) {
					double qq=0,*rp=rnoise+i+j*(m+2*k);
					for(int j2=-k;j2<=k;j2++) for(int i2=-k;i2<=k;i2++)
						qq+=rp[i2+j2*(m+2*k)]*exp(-(i2*i2+j2*j2)*fac);
					f.chi*=1+scale*qq;
				}
			}

			// Set staggered fields
			if(cc[ije]<4) {
				set_irwin_velocity(f,xs,ys);
				f.X=xs;f.Y=ys;
			}
		}
	}

	// Remove the random Gaussian field if it was used
	if(random_chi0) delete [] rnoise;

	fs->extrapolate();
}

void sbfrac::set_boundaries() {
	c_field *fp,*fm=fs->fm;
	double dx=fs->dx;
	int i;

	// Set bottom boundary
	for(i=1,fp=fm+1;i<m;i++,fp++) {
		set_irwin_velocity(*fp,-sim_size+dx*i,-sim_size);
		fp->X=2*fp[me].X-fp[2*me].X;
		fp->Y=2*fp[me].Y-fp[2*me].Y;
	}

	// Set top boundary
	for(i=1,fp=fm+me*m+1;i<m;i++,fp++) {
		set_irwin_velocity(*fp,-sim_size+dx*i,sim_size);
		fp->X=2*fp[-me].X-fp[-2*me].X;
		fp->Y=2*fp[-me].Y-fp[-2*me].Y;
	}

	// Set side boundaries
	for(i=0,fp=fm;i<me;i++) {
		set_irwin_velocity(*fp,-sim_size,-sim_size+dx*i);
		fp->X=2*fp[1].X-fp[2].X;
		fp->Y=2*fp[1].Y-fp[2].Y;fp+=m;
		set_irwin_velocity(*fp,sim_size,-sim_size+dx*i);
		fp->X=2*fp[-1].X-fp[-2].X;
		fp->Y=2*fp[-1].Y-fp[-2].Y;fp++;
	}
}

void sbfrac::set_irwin_velocity(c_field &f,double xs,double ys) {
	double r=sqrt(xs*xs+ys*ys),theta=r<1e-20?0:arg(xs,ys);
	f.u=Factor*sqrt(r)*(mix_i*((2*Kappa-1)*cos(theta*0.5)-cos(theta*1.5))
			  +mix_ii*((2*Kappa+3)*sin(theta*0.5)+sin(theta*1.5)));
	f.v=Factor*sqrt(r)*(mix_i*((2*Kappa+1)*sin(theta*0.5)-sin(theta*1.5))
			  +mix_ii*(-(2*Kappa-3)*cos(theta*0.5)-cos(theta*1.5)));
}

void sbfrac::random_chi_init(double *&rnoise,int &k,double &fac,double &scale) {
	int i,j;
	double qq=0;
	fac=random_chi0_len*fs->xsp;
	fac=0.5/(fac*fac);

	// Calculate how many gridpoints need to be scanned
	k=find_gaussian_cutoff(fac);
	printf("%d\n",k);
	fac*=2;
	for(i=0;i<k;i++) for(j=0;j<k;j++) qq+=(i==0?(j==0?1:2):(j==0?2:4))*exp(-(i*i+j*j)*fac);
	fac*=0.5;
	scale=random_chi0_pct*0.01/sqrt(qq);

	// Initialize the extended grid to be full of Gaussian noise
	int mee=m+2*k,meesq=mee*mee;
	double *rp=rnoise=new double[meesq];
	for(;rp<rnoise+(meesq-1);rp+=2) box_muller(*rp,rp[1]);
	if(rp!=rnoise+meesq) box_muller(*rp,qq);
}

inline int sbfrac::find_gaussian_cutoff(double fac) {
	const int max_sgrid=100;
	const double cutoff=1e-8;
	for(int k=1;k<max_sgrid;k++) if(exp(-k*k*fac)<cutoff) return k;
	fatal_error("Gaussian sampling too large",1);
	return 0;
}
