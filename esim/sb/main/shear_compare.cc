#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;

// The Bolztmann constant
const double kB=1.3806503e-23;

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc!=2) {
		fputs("Usage: ./shear_compare <case>\n",stderr);
		return 1;
	}

	// Check that the run number is in range
	int run=atoi(argv[1]);
	if(run<0||run>8) {
		fputs("Case number out of range\n",stderr);
		return 1;
	}

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

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Scaling factor from physical values
	const double sctab[9]={1e4,1e4,1e4,5e3,2.5e3,1.25e3,5e3,2.5e3,1.25e3};
	double sca=sctab[run];

	// Other parameters. Note that parameters are scaled to be larger than
	// realistic values, but allow for a direct--quasistatic comparison.
	const double visc=0.02;			// Viscous damping
	double t_scale=4.04645058986732e-06*sca;// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	double u_bdry=1e-7*sca;			// Wall speed (*)

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the two simulations
	const int ss[9]={160,160,160,160,160,160,160,160,160};
	shear_sim dsim(4*ss[run],ss[run],-4,4,-1,1,mu,K,visc,0,0,t_scale,tmult,adapt_fac,u_bdry,&stz,"");
	shear_sim qsim(4*ss[run],ss[run],-4,4,-1,1,mu,K,visc,0,0,t_scale,tmult,adapt_fac,u_bdry,&stz,"");

	// Set up fields
	const double chitab[9]={600,630,660,600,600,600,630,630,630};
	const double dchitab[9]={200,170,140,200,200,200,170,170,170};
	dsim.init_fields(1,chitab[run],dchitab[run]);
	qsim.init_fields(1,chitab[run],dchitab[run]);
	dsim.time=qsim.time=0;

	// Advance the two simulations forward in time
	double l2[4],target_time,time_interval=1000/sca;
	const double dt=dsim.dx*dsim.dx*dsim.tmult;
	const double nqtab[9]={5,5,5,5,5,5,80,80,80};
	double qdt=time_interval/nqtab[run];

	// Open the output file
	char buf[16];
	sprintf(buf,"sc%d.fe",run);
	FILE *fp=fopen(buf,"w");
	if(fp==NULL) {
		fputs("Can't open output file\n",stderr);
		return 1;
	}

	// Do an initial comparison of the fields
	dsim.l2_comparison(qsim,l2);
	fprintf(fp,"%g %.12g %.12g %.12g %.12g\n",dsim.time,*l2,l2[1],l2[2],l2[3]);

	for(int i=1;i<=2000;i++) {
		fflush(fp);

		// Direct simulation update
		target_time=time_interval*i;
		while(dsim.time+dt*(1+1e-11)<target_time) dsim.step_forward(dt);
		dsim.step_forward(target_time-dsim.time);

		// Quasi-static simulation update
		for(int j=0;j<nqtab[run];j++) qsim.step_forward_quasistatic(qdt);

		// Compare fields
		dsim.l2_comparison(qsim,l2);
		fprintf(fp,"%g %.12g %.12g %.12g %.12g\n",dsim.time,*l2,l2[1],l2[2],l2[3]);
	}
	fclose(fp);
}
