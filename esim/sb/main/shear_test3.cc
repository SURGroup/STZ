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

// This function returns a random double between 0 and 1
inline double rnd() {return double(rand())/RAND_MAX;}

void syntax_message() {
	fputs("Syntax: ./shear_test3 direct    (for direct simulation)\n"
	      "        ./shear_test3 qs        (for quasi-static simulation)\n",stderr);
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
	const double Delta=3500/TZ;
	const double Omega=133;
	const double eps0=0.3;
	const double chi_inf=950/TZ;
	const double theta=600/TZ;
	const double Omegaeps0=Omega*1e-30*eps0*1.1e9/(TZ*kB);

	// Output filenames
	const char dfn[]="sim_drnd",qfn[]="sim_qrnd";

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double visc=0.1;			// Viscous damping
	const double t_scale=0.2;		// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	const double u_bdry=0.001;		// Wall speed (*)

	// Make the output directory if it doesn't already exist
	mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	shear_sim sim(500,500,-2,2,-2,2,mu,K,visc,0,0,t_scale,tmult,adapt_fac,u_bdry,&stz,qs?qfn:dfn);
	sim.init_fields(0);

	// Code to read in a text file
	const int mm=50,nn=50;
	double *arr=new double[mm*nn],*arp=arr,*are=arr+mm*nn;

	// Read each value into the array
	while(arp<are) *(arp++)=rnd();

	// Initialize the chi field as the bicubic interpolation of the array
	sim.initialize_chi_bicubic(mm,nn,arr,580/TZ,40/TZ);
	delete [] arr;

	// Carry out the simulation using the selected simulation method
	qs?sim.solve_quasistatic(400,200,20):sim.solve(200,200);
}
