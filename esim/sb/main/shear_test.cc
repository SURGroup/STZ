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

void syntax_message() {
	fputs("Syntax: ./shear_test direct    (for direct simulation)\n"
	      "        ./shear_test qs        (for quasi-static simulation)\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check command-line arguments
	if(argc!=2) syntax_message();
	bool qs=true;
	if(strcmp(argv[1],"direct")==0) qs=false;
	else if (strcmp(argv[1],"qs")!=0) syntax_message();

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
	const char dfn[]="sim5d",qfn[]="sim5q";

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double sca=2e4;
	const double visc=0.02;			// Viscous damping
	const double t_scale=4.04645058986732e-06*sca; // Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	const double u_bdry=1e-7*sca;		// Wall speed (*)

	// Make the output directory if it doesn't already exist
	mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	shear_sim sim(480,160,-3,3,-1,1,mu,K,visc,0,0,t_scale,tmult,adapt_fac,u_bdry,&stz,qs?qfn:dfn);
	sim.init_fields(0,600,200);

	// Carry out the simulation using the selected simulation method
	qs?sim.solve_quasistatic(2e6/sca,400,10):sim.solve(2e6/sca,400);
}
