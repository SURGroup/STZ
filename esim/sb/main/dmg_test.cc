#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "dmg_sim.hh"

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
	const double a1=0;
	const double a2=0.6;
	const double E=70e9;			// Young's modulus for individual grains (Pa)
	const double gamma_g=1;			// Surface energy (J/m)
	const double C=1;			// O(1) constant
	const double V0=0.7;			// Initial volume fraction

	// Output filenames
	const char dfn[]="sdmg_d",qfn[]="sdmg_q";

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double ss=30e6*a2;		// Nominal stress scale (Pa)
	const double mu=110e6/ss;		// Rescaled shear modulus
	const double nu=0.3;			// Poisson ratio
	const double K=2*mu*(1+nu)/(3*(1-2*nu));// Rescaled bulk modulus
	const double rho0=1.6e3;		// Density (kg/m^3)
	const double cs=sqrt(mu*ss/rho0);	// Shear wave speed (m/s)

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double lscale=1;			// Length scale (m)
	const double sca=1;			// Artificial scaling
	const double visc=0.04;			// Viscous damping
	const double t_scale=lscale/cs*sca;	// Plasticity timescale (s) (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	const double u_bdry=0.2/cs*sca;		// Rescaled wall speed (*)
	const double aa_init=1e-4/lscale;	// Rescaled initial grain size
	printf("%g %g %g\n",u_bdry,t_scale,cs);

	// Make the output directory if it doesn't already exist
	mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_damage stz(a1,a2,mu,0.1,E/ss,gamma_g/(ss*lscale),V0,C,aa_init);

	// Initialize the simulation
	dmg_sim sim(1000,100,-0.5,0.5,-0.05,0.05,mu,K,visc,t_scale,tmult,adapt_fac,u_bdry,&stz,qs?qfn:dfn);
	sim.init_fields(0,0.03,0.07,30e6/ss,aa_init);

	// Carry out the simulation using the selected simulation method
	qs?sim.solve_quasistatic(0,20,200,4):sim.solve(0,2e6/sca,200);
}
