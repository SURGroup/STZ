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
const double TZ=21000;

// The Bolztmann constant
const double kB=1.3806503e-23;

int main() {

	// STZ model parameters (based on Vitreloy 1 BMG)
	const double c0=0.05;
	const double tau0=1e-13;
	const double kappa=4.2;
	const double Delta=3500/TZ;
	const double Omega=133;
	const double eps0=0.3;
	const double chi_inf=620/TZ; //950/TZ;
	const double theta=600/TZ;
	const double Omegaeps0=Omega*1e-30*eps0*1.1e9/(TZ*kB);

	// Output filenames
	const char fn[]="sim3t";

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

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double visc=0.02;			// Viscous damping
	const double chi_diff=0.;		// Effective temperature diffusion
	const double chi_len=0.;		// Dpl-mediated diffusion
	const double t_scale=0.2;		//0.2 Plasticity timescale (*)
	const double tmult=0.5;	//0.5		// Direct sim. timestep multiplier
	const double adapt_fac=0.01;		// Factor controlling plastic adaptivity
	const double u_bdry=0;			// 0.001 Wall speed (*)
	const double fx=0.5;			// Horizontal body force per unit volume

	// Make the output directory if it doesn't already exist
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_adam_original stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0,mu);

	// Initialize the simulation
	shear_sim sim(320,320,-2,2,-2,2,mu,K,visc,chi_diff,chi_len,t_scale,tmult,adapt_fac,u_bdry,fx,&stz,fflags,fn);

	// Code to read in a text file as the initial chi field
	ifstream file;
	file.open("initial_chi_pinaki");
	if(!file.is_open()) {
		cerr << "Error opening file." << endl;
		return -1;
	}

	// Read in the array dimensions
	int numRows, numCols;
	file >> numRows >> numCols;
	double *arr=new double[numRows*numCols];

	// Read each value into the array
	for(int i=0;i<numRows;i++) {
		for(int j=0;j<numCols; j++) {
			file >> arr[i+numRows*j];
		}
	}
	file.close();

	// Initialize the chi field as the bicubic interpolation of the array
	sim.initialize_chi_bicubic(numRows,numCols,arr,0.029,0.0095);
	delete [] arr;

	// Ramp up the applied force using the quasi-static solver
	sim.solve_quasistatic(100,100,10);

	// Deal with the plastic deformation using the direct solver
	sim.solve(1,100);
}
