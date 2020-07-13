#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

#include "common.hh"
#include "extra.hh"
#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=1;//21000;
const double TZ_=21000;                  // STZ formation energy, eZ/kB [K]
const double kB_=1.3806503e-23;           // Boltzmann constant [J/K]
const double kB_metal=8.617330350e-5;    // Boltzmann constant [eV/K]
const double eV_per_J=6.241509e18;       // Converts eV to Joules
const double Ez=TZ*kB_*eV_per_J;          // STZ formation Energy [eV]

// The Bolztmann constant
const double kB=1;

void syntax_message() {
    fputs("Syntax: ./shear_energy_Adam <run_type> <input_file> <beta> <u0> [additional arguments]\n\n"
          "<run_type> is either:\n"
           "     \"direct\" for direct simulation,or\n"
           "     \"qs\" for quasi-statc simulation\n\n"
           "<input_file> expects the relative path for a text file containing potential energy data\n\n"
           "<beta> is the initial energy scaling parameter [1/eV]\n\n"
           "<u0> is an energy offset value [eV]\n",stderr);
     exit(1);
}

int main(int argc,char **argv) {

	// Check command-line arguments
	if(argc<5) syntax_message();
	bool qs=true;
	if(strcmp(argv[1],"direct")==0) qs=false;
	else if(strcmp(argv[1],"qs") != 0) syntax_message();
    double beta=atof(argv[3]),u0=atof(argv[4]);

    printf("Beta is set to %g 1/eV\n", beta);
    printf("u0 is set to %g eV\n", u0);

	// STZ model parameters
	double c0=0.3;              //Work fraction 0.3
	double tau0=1e-13;			// Timescale (s)
	double chi_inf=2730;    //0.13 initially 2730
    double ep=10;           // STZ size    10
    double l_chi=4.01;      // 4.01
    
    char er_default[]="pe.Cu.MD.100.txt";
    char* endref=er_default;

    // Read additional command-line arguments to override default parameter values
    int i=5;
    while(i<argc) {
        if(read_arg(argv,"chi_inf",i,argc,chi_inf,"steady-state effective temperature"," K")) {}
        else if(read_arg(argv,"chi_len",i,argc,l_chi,"chi diffusion length scale"," Angstroms")) {}  
        else if(read_arg(argv,"c0",i,argc,c0,"plastic work fraction","")) {}
        else if(read_arg(argv,"ep",i,argc,ep,"STZ size","")) {}

        else if(se(argv[i],"endref")) {
            if(++i==argc) {
               fputs("Error reading command-line arguments\n",stderr);
               return 1;
            }
            printf("Reading final state from file %s\n",endref=argv[i]);
        } else {
            fprintf(stderr,"Command-line argument '%s' not recognized\n",argv[i]);
        }
        i++;
    }

    // If chi_inf was specified on the command line, then use that value.
    if(chi_inf>0) {
        printf("\nThe upper-limiting effective temperature was set using a user-defined value.\n"
               "The effective temperature is %g K.\n",chi_inf);
        chi_inf/=TZ_;
    } else {

        // Open final MD PE reference file and determine maximum PE
        FILE *fp=safe_fopen(endref,"r");
        double PE,maxPE=0;
        while(fscanf(fp,"%lf",&PE)) {
            if(PE>maxPE) maxPE=PE;
        }
        fclose(fp);

        // Compute value of chi_infinity from maximum potential energy
        chi_inf=0.95*beta*(maxPE-u0)*Ez/kB_metal/TZ;
        printf("\nThe maximum PE from the final snapshot is %g eV.\n"
               "For beta=%g and E0=%g the corresponding dimensionless value of"
               "chi_inf is %g. This is %g K.\n",maxPE,beta,u0,chi_inf,chi_inf*TZ);
    }

    // Output filenames
	const char dfn[]="sct_d.out",qfn[]="sct_q.out";

    // Flags governing the files to output
    const unsigned int fflags=4|8|16|32|64|128|256|8192|32768;
    //4|8|16|32|64|128|256|8192|32678;

	// Elasticity related parameters
	const double mu_phys=20;                // Shear modulus (GPa)
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
    printf("Bulk modulus is %g\n",K);

	// Make the output directory if it doesn't already exist
	mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Set up plasticity model
	//const double ep=10;
	stz_dynamics_adam stz_a(TZ,chi_inf,tau0,ep,c0);

	// Initialize the simulation
	shear_sim sim(32,32,0,1,0,1,mu,K,visc,chi_diff,chi_len,t_scale,tmult,adapt_fac,u_bdry,0,&stz_a,fflags,qs?qfn:dfn);
    sim.init_fields(0,580,220);
	// sim.initialize_tracers(32,32);
    	
	// Code to read in a text file as the initial chi field
	// ifstream file;
	// file.open("sig_50_rand2");//strain_diff_8x8
	// if(!file.is_open()) {
	//	cerr << "Error opening file." << endl;
	//	return -1;
	// }

	// Read in the array dimensions
	// int numRows, numCols;
	// file >> numRows >> numCols;
	// double *arr=new double[numRows*numCols];

	// Read each value into the array
	// for(int i=0;i<numRows;i++){ for(int j=0;j<numCols; j++){
    //		file >> arr[j+numRows*i];
	//       }
	// }
	// file.close();

	// Initialize the chi field as the bicubic interpolation of the array
	// sim.initialize_chi_bicubic(numRows,numCols,arr,0,1);
	// delete [] arr;
    // sim.initialize_random(2,570,10);
	
    // Open the input file and read in the grid dimensions and scaling
    read_chi_from_file(argv[2],sim,u0,beta*Ez/kB_metal,TZ); 

    // Carry out the simulation using the selected simulation method
	//solve_quasistatic(double duration,int frames,int steps)
	//qs?sim.solve_quasistatic(drn,100,40):sim.solve(drn,100);
    qs?sim.solve_quasistatic(drn,100,60):sim.solve(drn,100);
}
