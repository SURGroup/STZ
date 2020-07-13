#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=1;//21000;

// The Bolztmann constant
const double kB=1;

void syntax_message() {
	fputs("Syntax: ./shear_test4 direct    (for direct simulation)\n"
	      "        ./shear_test4 qs        (for quasi-static simulation\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check command-line arguments
	if(argc!=2) syntax_message();
	bool qs;
	if(strcmp(argv[1],"direct")==0) qs=false;
	else if (strcmp(argv[1],"qs")==0) qs=true;
	else syntax_message();

	// STZ model parameters
	const double c0=1;
	const double tau0=1e-13;			// Timescale (s)
	const double chi_inf=0.15;

	// Output filenames
	const char dfn[]="sim4d",qfn[]="sim4q";

	// Flags governing the files to output
	// 1 - u	256 - dev
	// 2 - v	512 - X
	// 4 - p	1024 - Y
	// 8 - q	2048 - QS diagnostics
	// 16 - s	4096 - Dplastic
	// 32 - tau	8192 - Lagrangian fields
	// 64 - chi	16384 - Dtotxy
	// 128 - tem
	const unsigned int fflags=64|128|256;

	// Elasticity related parameters
	const double mu_phys=25;                // Shear modulus (GPa)
	const double nu=0.35;                   // Poisson ratio
	const double s_y=0.85;                  // Yield stress (GPa)
	const double mu=mu_phys/s_y;            // Shear modulus
	const double K=2*mu*(1+nu)/(3*(1-2*nu));// Bulk modulus

	// Other parameters
	const double l_scale=4.01e-8;           // Length scale (m)
	const double rho0=6.125e3;		// Density (kg/m^3)
	const double c_s=sqrt(mu_phys*1e9/rho0);// Shear wave speed (m s^{-1})
	//const double c_s=sqrt(mu_phys*1e22/rho0);// Shear wave speed (m s^{-1})
	const double t_scale=l_scale/c_s;	// Time scale (s)
	const double adapt_fac=0.002;		// Factor controlling plastic adaptivity
	const double visc=0;			// Viscous damping
	const double cd_phys=0;//4.01e-10;		// Chi diffusion (m^2 s^{-1})
	const double chi_diff=cd_phys*t_scale/(l_scale*l_scale); // Chi diffusion
	const double chi_len=0.01;		// Dpl-mediated diffusion
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double s_rate=1e8;		// Strain rate (s^{-1})
	const double u_bdry_phys=0.5*s_rate*l_scale; // Wall speed (m/s)
	const double u_bdry=u_bdry_phys/c_s;	// Wall speed
	const double drn_phys=0.5e-8;		// Simulation duration (s)
	const double drn=drn_phys/t_scale;	// Simulation duration

	printf("%g %g %g %g %g\n",c_s,t_scale,chi_diff,u_bdry,drn);

	// Make the output directory if it doesn't already exist
	mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Set up plasticity model
	const double ep=10;
	stz_dynamics_adam stz_a(TZ,chi_inf,tau0,ep,c0);

	// Initialize the simulation
	shear_sim sim(200,200,0,1,0,1,mu,K,visc,chi_diff,chi_len,t_scale,tmult,adapt_fac,u_bdry,0,&stz_a,fflags,qs?qfn:dfn);

	// Code to read in a text file as the initial chi field
	ifstream file;
	file.open("window_50x50");//strain_diff_8x8
	if(!file.is_open()) {
		cerr << "Error opening file." << endl;
		return -1;
	}

	// Read in the array dimensions
	int numRows, numCols;
	file >> numRows >> numCols;
	double *arr=new double[numRows*numCols];

	// Read each value into the array
	for(int i=0;i<numRows;i++){ for(int j=0;j<numCols; j++){
		file >> arr[j+numRows*i];
	        }
	}
	file.close();

	// Initialize the chi field as the bicubic interpolation of the array
	sim.initialize_chi_bicubic(numRows,numCols,arr,0,1);
	//delete [] arr;

	// Carry out the simulation using the selected simulation method
	//solve_quasistatic(double duration,int frames,int steps)
	qs?sim.solve_quasistatic(drn,100,40):sim.solve(drn,100);

}
