#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "shear_sim.hh"

#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
inline double wtime() {return 0;}
#endif

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;

// The Bolztmann constant
const double kB=1.3806503e-23;

void output_tractions(shear_sim &sim,FILE *fp,const double t_scale,const double u_bdry) {
	double bx,by,tx,ty,Q,max_qs,&t=sim.time,bpos=sim.bdry_pos();
	sim.tractions(bx,by,tx,ty);
	sim.qs_measure(Q,max_qs);
	fprintf(fp,"%g %g %g  %.12g %.12g %.12g %.12g  %.12g %.12g\n",t,t*t_scale,bpos,bx,by,tx,ty,Q,max_qs);
	fflush(fp);
	//printf("T. info.: %g %g %g  %g %g %g %g  %g %g\n",t,t*t_scale,bpos,bx,by,tx,ty,Q,max_qs);
}

int main(int argc,char **argv) {

	if(argc!=4) {
		fputs("Syntax: ./shear_random <smooth scale> <initial chi> <chi variation>\n",stderr);
		return 1;
	}
	double smsca=atof(argv[1]),ichi=atof(argv[2]),vchi=atof(argv[3]);

	// STZ model parameters (based on Vitreloy 1 BMG)
	const double c0=0.4;
	const double tau0=1e-13;
	const double kappa=4.2;
	const double Delta=8000/TZ;
	const double Omega=300;
	const double eps0=0.3;
	const double chi_inf=900/TZ;
	const double theta=400/TZ;
	const double Omegaeps0=Omega*1e-30*eps0*1.1e9/(TZ*kB);

	// Output filenames
	char fn[256];
	sprintf(fn,"simh_%s_%sK_%sK",argv[1],argv[2],argv[3]);

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=90e9;			// Young's modulus
	const double nu=0.375;			// Poisson ratio
	const double s_y=1.75e9;		// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus
	const double den=7.7e3;			// Density (kg/m^3)
	const double l_scale=0.001;		// Length scale

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double visc=0.04;			// Viscous damping
	const double t_scale=l_scale*sqrt((2*(1+nu)*den)/E);// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	const double s_rate=1e-4;		// Strain rate (1/s)
	const double u_bdry=s_rate*t_scale;	// Boundary velocity
	const double max_strain=0.3;		// Maximum strain
	const double chi_diff=0.;		// Chi diffusion
	const double chi_len=0.01;		// Dpl-mediated diffusion

	// Make the output directory if it doesn't already exist
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	shear_sim sim(2000,1000,-2,2,-1,1,mu,K,visc,chi_diff,chi_len,t_scale,tmult,adapt_fac,u_bdry,&stz,fn);
	sim.initialize_random(smsca,ichi,vchi);

	// Carry out the simulation using the selected simulation method
	sim.time=0;
	int steps=40;
	const int frames=500;
	const double t_end=max_strain/u_bdry,time_interval=t_end/frames;
	double t0,t1,t2;

	// Open special diagnostic file
	char fn2[260];
	sprintf(fn2,"%s.trf",fn);
	FILE *fp=safe_fopen(fn2,"w");
	output_tractions(sim,fp,t_scale,u_bdry);

	// Output the initial fields
	sim.write_files(0);
	puts("# Frame 0");
	t0=wtime();

	for(int k=1;k<=frames;k++) {

		for(int j=0;j<steps;j++) {
			sim.step_forward_quasistatic(time_interval/steps);
			output_tractions(sim,fp,t_scale,u_bdry);
		}
		t1=wtime();

		// Output the fields
		sim.write_files(k);

		// Print diagnostic information
		t2=wtime();
		printf("# Frame %d [%d, %.8g s, %.8g s] {MG %.2f}\n",k,steps,t1-t0,t2-t1,sim.avg_iters());
		t0=t2;
	}

	fclose(fp);
}
