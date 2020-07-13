#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;

// The Bolztmann constant
const double kB=1.3806503e-23;

void syntax_message() {
    fputs("Syntax: ./shear_data <type> <input_file>\n\n"
          "Type is \"direct\" for direct simulation\n"
          "        \"qs\" for quasi-static simulation\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

    // Check command-line arguments
    if (argc != 3) syntax_message();
    bool qs = true;
    if (strcmp(argv[1], "direct") == 0) qs = false;
    else if (strcmp(argv[1], "qs") != 0) syntax_message();

    // STZ model parameters (based on Vitreloy 1 BMG)
    const double s_y = 0.85e9;              // Yield stress (Pa)
    const double tau0 = 1e-13;

    // Output filenames
    const char dfn[] = "sct_d.out", qfn[] = "sct_q.out";

    // Output fields and miscellaneous flags. 1-u, 2-v, 4-p, 8-q, 16-s, 32-tau,
    // 64-chi, 128-tem, 256-dev, 512-X, 1024-Y, 2048-(total strain components),
    // 4096-(total strain invariants), 8192-(Lagrangian tracer output).
    // const unsigned int fflags=128|256|512|1024|8192;

    const unsigned int fflags=4|8|16|32|128|2048;

    // Elasticity related parameters (based on Vitreloy 1 BMG)
    const double E = 80e9;                 // Young's modulus (Pa)
    const double nu = 0.35;                 // Poisson ratio
    //const double K = E/(3*(1-2*nu)*s_y);    // Rescaled bulk modulus
    //const double mu = E/(2*(1+nu)*s_y);     // Rescaled shear modulus
    const double K = 143.23;
    const double mu = 20;

    // Other parameters. Note that parameters labeled (*) are far larger
    // than realistic values, but allow for a direct--quasistatic
    // comparison.
    const double le = 1e-2;                 // Length scale (m)
    const double rho0 = 6125;               // Density (kg/m^3)
    const double sca = 2e4;
    const double t_scale = le/sqrt(mu*s_y/rho0)*sca; // Plasticity timescale (*)
    const double visc = 0.02;               // Viscous damping
    // const double chi_len=0.06;              // Dpl-mediated diffusion
    const double chi_len=4.01;
    const double adapt_fac = 1e-4;          // Adaptivity factor(initial values 1e-4)

    // Set parameters for either a periodic or non-periodic test
    const bool y_prd=true;
    double u_bdry, lamb, tf;
    if(y_prd) {

        // Periodic test
        u_bdry=0.;
        lamb=1e-3;tf=0.5/lamb;
    } else {

        // Non-periodic test
        u_bdry=1e-7*sca;
        lamb=0;tf=3e6/sca;
    }

    // Make the output directory if it doesn't already exist
    mkdir(qs? qfn : dfn, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Initialize an STZ plasticity model
    // The Bouchbinder-Rycroft model from 2012 PRL, 2015 JCP, 2016 PRA papers.
    /*const double c0 = 0.4;
    const double kappa = 4.2;
    const double Delta = 8000/TZ;
    const double Omega = 300.*22/17;        // Fix typo in 2015 JCP
    const double eps0 = 0.3;
    const double chi_inf = 900/TZ;
    const double theta = 400/TZ;
    const double Omegaeps0 = Omega*1e-30*eps0*s_y/(TZ*kB);
    stz_dynamics_linear_athermal stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);*/

    // The model from Hinkle et al. (PRE, 2017)
	const double chi_inf=0.15; 
    const double ep=10;
    const double c0=1;
    stz_dynamics_adam stz(TZ,chi_inf,tau0,ep,c0);

    // Initialize the simulation
    shear_sim sim(32, 32, 0, 400, 0, 400, mu, K, visc, chi_len, t_scale, adapt_fac, u_bdry, lamb, &stz, y_prd, fflags, qs? qfn : dfn);
    sim.init_fields(0, 580, 220);
    //sim.initialize_tracers(64,64);

    // Open the input file and read in the grid dimensions and scaling
 /*   FILE *fp=safe_fopen(argv[2],"r");
    int mm,nn;
    double chi_base,chi_scale;
    if(fscanf(fp,"%d %d %lg %lg",&mm,&nn,&chi_base,&chi_scale)!=4) {
        fputs("Error reading header information\n",stderr);
        return 1;
    }

    // Check that the dimensions make sense, and allocate memory
    if(mm<=0 || nn<=0) {
        fputs("Grid dimensions are invalid\n",stderr);
        return 1;
    }
    double *f=new double[mm*nn];

    // Read in the chi values from the file
    for(int j=0;j<nn;j++) {
        for(int i=0;i<mm;i++) {
            if(fscanf(fp,"%lg",f+i+mm*j)!=1) {
                fputs("Error reading chi information\n",stderr);
                return 1;
            }
        }
    }
    printf("# chi grid dimensions : %d by %d\n"
           "# chi base scale      : %g K\n"
           "# chi range scale     : %g K\n",mm,nn,chi_base,chi_scale);
    sim.initialize_chi_bicubic(mm,nn,f,chi_base/TZ,chi_scale/TZ);

    // Free the dynamically allocated memory
    delete [] f;*/

    // Initialize random chi field
    sim.initialize_random(2,570,10);

    // Carry out the simulation using the selected simulation method
    qs?sim.solve_quasistatic(tf,100,10):sim.solve(tf,160);
}
