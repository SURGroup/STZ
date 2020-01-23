#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "bd_sim.hh"
#include "extra_force.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;

// The Bolztmann constant
const double kB=1.3806503e-23;

void syntax_message() {
	fputs("Syntax: ./bd_test <case>\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check command-line arguments
	if(argc!=2) syntax_message();
	int run=atoi(argv[1]);
	if(run<1||run>5) {
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

	// Output filenames
	char fn[16];
	sprintf(fn,"bdry%d",run);
	const unsigned int fflags=4|128|256|512|1024|1048576;

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double sca=1;
	const double visc=0.02;			// Viscous damping
	const double t_scale=4.04645058986732e-06*sca;// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor

	// Make the output directory if it doesn't already exist
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	if(run==1) {
		bd_sim sim(512,256,-2,2,-1,1,mu,K,visc,0,0.001,t_scale,tmult,adapt_fac,&stz,fn,fflags|12288);
		extra_force1 eforce(sim);
		sim.init_bar(1);
		sim.solve_quasistatic(100*1e4,400,2);
	} else if(run==2) {
		bd_sim sim(512,1024,-0.5,0.5,-1,1,mu,K,visc,0,0.008,t_scale,tmult,adapt_fac,&stz,fn,fflags);
		extra_force2 eforce(sim);
		sim.init_bar(2);
		sim.solve_quasistatic(50*1e4,200,4);
		sim.solve(10,200);
	} else if(run==3) {
		bd_sim sim(512,256,-2,2,-1,1,mu,K,visc,0,0.001,t_scale,tmult,adapt_fac,&stz,fn,fflags|12288);
		extra_force3 eforce(sim);
		sim.init_bar(1);
		sim.solve_quasistatic(100*1e4,400,2);
	} else if(run==4) {
		bd_sim sim(150,150,-1.6,1.6,-1.6,1.6,mu,K,visc*10,0,0,t_scale,tmult,adapt_fac,&stz,fn,1<<21);
		sim.init_bar(3);
		sim.solve(200,10000);
	} else {
		bd_sim sim(150,150,-1.6,1.6,-1.6,1.6,mu,K,visc*10,0,0,t_scale,tmult,adapt_fac,&stz,fn,1<<21);
		sim.init_bar(3);
		sim.wrong_spin=true;
		sim.solve(200,10000);
	}
}
