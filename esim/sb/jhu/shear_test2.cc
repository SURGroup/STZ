#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;		// "STZ Formation Energy"

// The Bolztmann constant
const double kB=1.3806503e-23;

void syntax_message() {
	fputs("Syntax: ./shear_test2 direct    (for direct simulation)\n"
	      "        ./shear_test2 qs        (for quasi-static simulation\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check command-line arguments
	if(argc!=2) syntax_message();
	bool qs;
	if(strcmp(argv[1],"direct")==0) qs=false;
	else if (strcmp(argv[1],"qs")==0) qs=true;
	else syntax_message();

	// STZ model parameters (based on Vitreloy 1 BMG)
	const double c0=0.4;			// Fraction of plastic work that goes to increasing chi (set to 0.4 in Rycroft paper)
	const double tau0=1e-13;		// Sets the rate at which STZs 'flip', "molecular vibration timescale" (generally set to 1e-13)
	const double kappa=4.2;			// Additional STZ model parameter
	const double Delta=8000/TZ; 		// Nondimensionalized "Typical Activation Barrier"
	const double Omega=300;			// Activation Volume
	const double eps0=0.3;			// Typical local strain at transition from elastic to plastic deformation
	const double chi_inf=900/TZ;		// Nondimensionalized saturated effective temperature
	const double theta=400/TZ;
	const double Omegaeps0=Omega*1e-30*eps0*1.1e9/(TZ*kB);

	// Output filenames
	const char dfn[]="sim2d.out",qfn[]="sim2q.out";

	// Flags governing the files to output
	// 1 - u	256 - dev
	// 2 - v	512 - X
	// 4 - p	1024 - Y
	// 8 - q	2048 - QS diagnostics
	// 16 - s	4096 - Dplastic
	// 32 - tau	8192 - Lagrangian fields
	// 64 - chi	16384 - Dtotxy
	// 128 - tem
	const unsigned int fflags=128|256|512|1024|4096;

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double sca=1e4;			// Scaling parameter from Chris's input
	const double visc=0.0001;		// Viscous damping
	const double chi_diff=0.;		// Effective temperature diffusion
	const double chi_len=0.06;		// Dpl-mediated diffusion
	const double t_scale=4.04645058986732e-06*sca;		// 0.2 Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Factor controlling plastic adaptivity
	const double u_bdry=1e-7*sca;		// 0.001 Wall speed (*)
	const double fx=0.0;			// Horizontal body force per unit volume

	// Make the output directory if it doesn't already exist
	mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	shear_sim sim(200,200,-1,1,-1,1,mu,K,visc,chi_diff,chi_len,t_scale,tmult,adapt_fac,u_bdry,fx,&stz,fflags,qs?qfn:dfn);

	// Read random field parameters from file
	double rf_mean, rf_std, rf_len;
	FILE *fp=safe_fopen("random_field_input","r");
	if(fscanf(fp,"%lg %lg %lg",&rf_mean,&rf_std,&rf_len)!=3) {
		fputs("Error reading parameters from file\n",stderr);
		return 1;
	}
	fclose(fp);

	// Initialize the effective temperature field
	sim.init_fields(0,rf_mean,rf_std);
	sim.initialize_tracers();
	sim.initialize_random(rf_len,rf_mean,rf_std);

	// Carry out the simulation using the selected simulation method
	qs?sim.solve_quasistatic(200,24,20):sim.solve(200,24);
}
