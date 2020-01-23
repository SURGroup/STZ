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
	fputs("Syntax: ./bd_test\n",stderr);
	exit(1);
}

int main() {

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
	const char fn[16]="bp.out";
	const unsigned int fflags=3|4|128|256|512|1024|1048576;

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double sca=1e4;			// Scaling factor from physically realistic value
	const double visc=0.02;			// Viscous damping
	const double chi_len=0.05;		// The chi diffusion length scale
	const double t_scale=4.04645058986732e-06*sca;// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	const double wallx=2;			// Walls initial boundary (originally 1)
	const double wallu=1e-7*sca;		// Walls drawing velocity (originally 1e-7*sca)

	// Make the output directory if it doesn't already exist
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	//bd_sim sim(512,256,-2,2,-1,1,mu,K,visc,0.001,t_scale,tmult,adapt_fac,false,&stz,fn,fflags|12288);
	bd_sim sim(1024,256,-4,4,-1,1,mu,K,visc,chi_len,0.001,t_scale,tmult,adapt_fac,&stz,fn,fflags|12288);

	// Set wall position and velocity
	sim.wallx=wallx;
	sim.wallu=wallu;

	// Wall boundary condition - set to either "clamped", "co_thin", or "incompr"
	sim.wall_bc=co_thin;

	// Initialize the simulation fields and solve
	sim.init_bar();
	sim.solve_quasistatic(3e7/sca,3000,2);
}
