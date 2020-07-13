#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "extra.hh"
#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;                   // STZ formation energy, eZ/kB [K]
const double kB=1.3806503e-23;           // Boltzmann constant [J/K]
const double kB_metal=8.617330350e-5;    // Boltzmann constant [eV/K]
const double eV_per_J=6.241509e18;       // Converts eV to Joules
const double Ez=TZ*kB*eV_per_J;          // STZ formation Energy [eV]

void syntax_message() {
    fputs("Syntax: ./shear_no_MD > [additional arguments]\n\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

    // Check command-line arguments
    bool qs=true;
    
    // Set default values of parameters
    double rho0=7270.7;           // Density (kg/m^3)
    double s_y=0.85;              // Yield stress (GPa)
    double Delta_=2000;           // Activation temperature (K) 
    double Omega=140;              // Activation volume (A^3)
    double l_chi=1;
    double c0=0.4;                // Plastic work fraction
    double eps0=0.3;              // Local strain
    double mu_=23.9;              // Shear modulus (GPa)
    double K_=110.42;             // Bulk modulus (GPa)
    double theta_=100;            // Bath temperature (K)
    double chi_inf_=2730;         // Upper limiting effective temperature (K)
    double kappa_=0;
	double T_scale_=4;            // Effective temperature scaling factor (?)
	double T_mean_=400;           // Effective temperature mean value
	double T_std_=50; 	          // Effective temperature standard deviation

    // Read additional command-line arguments to override default parameter values
    int i=1;
    while(i<argc) {
        if(read_arg(argv,"chi_inf",i,argc,chi_inf_,"steady-state effective temperature"," K")) {}
        else if(read_arg(argv,"theta",i,argc,theta_,"bath temperature"," K")) {}
        else if(read_arg(argv,"sy",i,argc,s_y,"yield stress"," GPa")) {}
        else if(read_arg(argv,"Delta",i,argc,Delta_,"activation barrier"," K")) {}
        else if(read_arg(argv,"Omega",i,argc,Omega,"activation volume"," Angstroms^3")) {}
        else if(read_arg(argv,"chi_len",i,argc,l_chi,"chi diffusion length scale"," Angstroms")) {}
        else if(read_arg(argv,"c0",i,argc,c0,"plastic work fraction","")) {}
        else if(read_arg(argv,"mu",i,argc,mu_,"shear modulus"," GPa")) {}
        else if(read_arg(argv,"K",i,argc,K_,"bulk modulus"," GPa")) {}
        else if(read_arg(argv,"rho",i,argc,rho0,"density"," kg/m^3")) {}
        else if(read_arg(argv,"eps0",i,argc,eps0,"local strain","")) {}
        else if(read_arg(argv,"kappa",i,argc,kappa_,"diffusion prefactor","")) {}
        else if(read_arg(argv,"T_scale",i,argc,T_scale_,"Effective temperature scale","")) {}
        else if(read_arg(argv,"T_mean",i,argc,T_mean_,"Effective temperature mean value","")) {}
        else if(read_arg(argv,"T_std",i,argc,T_std_,"Effective temperature standard deviation","")) {}
         else {
            fprintf(stderr,"Command-line argument '%s' not recognized\n",argv[i]);
        }
        i++;
    }


    // Rescaled elasticity parameters
    const double K=K_/s_y;                    // Rescaled bulk modulus (--)
    const double mu=mu_/s_y;                  // Rescaled Shear modulus (--)

    // STZ model parameters (based on Vitreloy 1 BMG)
    const double tau0=1e-13;                  // Vibration timescale (s)
    const double kappa=kappa_*1e30;
    const double Delta=Delta_/TZ;             // Activation barrier
    const double theta=theta_/TZ;             // Bath Temperature, (--)
	const double chi_inf=chi_inf_/TZ;         //     
    const double Omegaeps0=Omega*1e-30*eps0*s_y*1e9/(TZ*kB);
	const double T_scale=T_scale_;
	const double T_mean=T_mean_;
	const double T_std=T_std_;

    // Output filenames
    const char dfn[]="sct_d.out",qfn[]="sct_q.out";
    //const char* dfn=dfn_;
    //const char* qfn=qfn_;

    // Output fields and miscellaneous flags. 1-u,2-v,4-p,8-q,16-s,32-tau,
    // 64-chi,128-tem,256-dev,512-X,1024-Y,2048-(total strain components),
    // 4096-(total strain invariants),8192-(Lagrangian tracer output).
    const unsigned int fflags=32|128|512|1024|2048|4096;

    // Other parameters. The scale factor applies a scaling to the rate of
    // plasticity and the applied strain. It is used to allow for comparison
    // between the explicit and quasi-static models.
    const double le=1e-9;                   // Length scale (m)
    const double sca=2e4;                   // Scale factor [only used in non-periodic test]
    const double visc=0.02;                 // Viscous damping
    // const double chi_len=l_chi*1e-10;       // Dpl-mediated diffusion (m)
    const double chi_len=l_chi;       // Dpl-mediated diffusion (m)
    const double adapt_fac=2e-3;

    // MD simulation-based continuum parameters
    const int x_grid=32;              // Grid points in x-direction
    const int y_grid=32;              // Grid points in y-direction
    const double x_beg=0.0;           // x-origin (m)
    const double x_end=4e-8;          // x-terminus (m)
    const double y_beg=0.0;           // y-origin (m)
    const double y_end=4e-8;          // y-terminus (m)
    const double gam_max=0.5;         // MD max strain

    // Compute conversion factor from simulation time to physical time. This
    // enters into the plasticity model. This value is determined automatically
    // from the length scale, density, and shear modulus in order to make shear
    // waves move at speed 1 in the simulation units.
    const double t_scale=le/sqrt(mu*s_y/rho0);
    printf("One simulation length unit is %g nm\n"
           "One simulation time unit is %g ns\n",le*1e9,t_scale*1e9);

    // Set parameters for either a periodic or non-periodic test
    const bool y_prd=true;
    double u_bdry,lamb,lamb_phys,tf;
    if(y_prd) {

        // Periodic test
        u_bdry=0.;
        lamb_phys=1e8;                // Strain rate (1/s)
        lamb=lamb_phys*t_scale;       // Strain rate (sim. units)
        tf=gam_max/lamb;              // Final time (sim. units)
    } else {

        // Non-periodic test
        u_bdry=1e-7*sca;
        lamb=0;tf=3e6/sca;
    }

    // Make the output directory if it doesn't already exist
    mkdir(qs?qfn:dfn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Initialize an STZ plasticity model
    //stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);
    stz_dynamics_linear stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,rho0);

    // Initialize the simulation
    shear_sim sim(x_grid,y_grid,x_beg/le,x_end/le,y_beg/le,y_end/le,mu,K,
                  visc,chi_len,t_scale,adapt_fac,u_bdry,lamb,&stz,y_prd,fflags,qs? qfn : dfn);
    //sim.init_fields(0,580,220);

    // Open the input file and read in the grid dimensions and scaling
    //
    //
    sim.initialize_random(T_scale,T_mean,T_std);

    // Carry out the simulation using the selected simulation method
    int n_frames=100,steps=10;
    printf("Simulation time unit   : %g s\n",t_scale);
    printf("Final time             : %g (%g s) \n",tf,tf*t_scale);
    printf("Quasi-static step size : %g \n",tf/(n_frames*steps));
    qs?sim.solve_quasistatic(tf,n_frames,steps):sim.solve(tf,160);
}
