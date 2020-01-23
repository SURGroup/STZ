#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000; // STZ formation energy, eZ/kB [K]
const double kB=1.3806503e-23; // Boltzman constant, [J/K]
const double kB_metal = 8.617330350e-5; // Boltzmann constant, [eV/K]
const double eV_per_J = 6.241509e18; // Converts eV to Joules
const double Ez = TZ * kB * eV_per_J; // STZ formation Energy, [eV]
// const double Ez = 2.0;

// double beta, u0;

void syntax_message() {
    fputs("Syntax: ./shear_energy <run_type> <input_file> <beta> <u0>\n\n"
          "<run_type> is either:\n"
          "     \"direct\" for direct simulation, or\n"
          "     \"qs\" for quasi-static simulation. \n\n"
          "<input_file> expects the relative path for a text file containing potential energy data. \n\n"
          "<beta> is the initial energy scaling parameter, [1/eV]. \n\n"
          "<u0> is an energy offset value, [eV].\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

    // Check command-line arguments
    if (argc != 5) syntax_message();
    bool qs = true;
    if (strcmp(argv[1], "direct") == 0) qs = false;
    else if (strcmp(argv[1], "qs") != 0) syntax_message();
    double beta, u0;
    beta = strtod(argv[3],NULL);
    u0 = strtod(argv[4],NULL);

    // STZ model parameters (based on Vitreloy 1 BMG)
    const double s_y = 0.85e9;              // Yield stress (Pa)
    const double c0 = 0.4;
    const double tau0 = 1e-13;
    //const double tau0 = 1e-4;
    const double kappa = 4.2;
    const double Delta = 8000/TZ;
    const double Omega = 300.*22/17;        // Fix typo in 2015 JCP
    const double eps0 = 0.3;
    //const double chi_inf = 0.13*Ez/kB_metal/TZ; // A.H. Dissertation, pp. 94
    const double chi_inf = 900/TZ;
    const double theta = 400/TZ; // Bath Temperature, [--]
    const double Omegaeps0 = Omega*1e-30*eps0*s_y/(TZ*kB);

    // Output filenames
    const char dfn[] = "sct_d.out", qfn[] = "sct_q.out";

    // Output fields and miscellaneous flags. 1-u, 2-v, 4-p, 8-q, 16-s, 32-tau,
    // 64-chi, 128-tem, 256-dev, 512-X, 1024-Y, 2048-(total strain components),
    // 4096-(total strain invariants), 8192-(Lagrangian tracer output).
    const unsigned int fflags=128|256|512|1024|4096|8192;

    // Elasticity related parameters (based on Vitreloy 1 BMG)
    const double E = 101e9;                 // Young's modulus (Pa)
    const double nu = 0.35;                 // Poisson ratio
    const double K = E/(3*(1-2*nu)*s_y);    // Rescaled bulk modulus
    const double mu = E/(2*(1+nu)*s_y);     // Rescaled shear modulus

    // Other parameters. Note that parameters labeled (*) are far larger
    // than realistic values, but allow for a direct--quasistatic
    // comparison.
    //const double le = 1e-2;                 // Length scale (m)
    const double le = 500;                 // Length scale (m)
    const double rho0 = 6125;               // Density (kg/m^3)
    const double sca = 2e4;
    const double t_scale = le/sqrt(mu*s_y/rho0)*sca; // Plasticity timescale (*)
    const double visc = 0.02;               // Viscous damping
    const double chi_len=0.06;              // Dpl-mediated diffusion
    const double adapt_fac = 1e-4;          // Adaptivity factor

    // MD Simulation-Based Continuum Parameters
    const int x_grid = 32;              // grid-points in x-direction
    const int y_grid = 32;              // grid-points in y-direction
    const double x_beg = 0.0;           // x-origin
    const double x_end = 400.0;         // x-terminus
    const double y_beg = 0.0;           // y-origin
    const double y_end = 400.0;         // y-terminus
    const double gam_dot = 1e-4;        // MD strain rate
    //const double gam_dot = 1e-16;        // MD strain rate
    const double gam_max = 0.5;         // MD max strain

    // Set parameters for either a periodic or non-periodic test
    const bool y_prd=true;
    double u_bdry, lamb, tf;
    if(y_prd) {

        // Periodic test
        u_bdry=0.;
        lamb=gam_dot;
        tf=gam_max/lamb;
    } else {

        // Non-periodic test
        u_bdry=1e-7*sca;
        lamb=0;tf=3e6/sca;
    }

    // Make the output directory if it doesn't already exist
    mkdir(qs? qfn : dfn, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Initialize an STZ plasticity model
    stz_dynamics_linear_athermal stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);

    // Initialize the simulation
    shear_sim sim(x_grid, y_grid, x_beg, x_end, y_beg, y_end, mu, K, visc, chi_len, t_scale, adapt_fac, u_bdry, lamb, &stz, y_prd, fflags, qs? qfn : dfn);
    sim.init_fields(0, 580, 220);

    // Open the input file and read in the grid dimensions and scaling
    FILE *fp=safe_fopen(argv[2],"r");
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

    // Read in the energy values from the file
    for(int j=0;j<nn;j++) {
        for(int i=0;i<mm;i++) {
            if(fscanf(fp,"%lg",f+i+mm*j)!=1) {
                fputs("Error reading energy information\n",stderr);
                return 1;
            }
        }
    }

    // Map energy values onto effective temperature, chi
    for(int j=0; j<nn; j++) {
        for(int i=0; i<mm; i++) {
            f[i+mm*j] = beta*(f[i+mm*j]-u0)*Ez/kB_metal;
        }
    }

    printf("# chi grid dimensions : %d by %d\n"
           "# chi base scale      : %g K\n"
           "# chi range scale     : %g K\n",mm,nn,chi_base,chi_scale);
    sim.initialize_chi_bicubic(mm,nn,f,chi_base/TZ,chi_scale/TZ);

    // Free the dynamically allocated memory
    delete [] f;

    // Carry out the simulation using the selected simulation method
    // qs?sim.solve_quasistatic(tf,n_frames,5):sim.solve(tf,160);
    //double duration = 0.5;
    double duration = tf;
    int n_frames = 100;
    int steps = 10;
    // qs?sim.solve_quasistatic(tf,n_frames,1):sim.solve(tf,160);
    qs?sim.solve_quasistatic(duration,n_frames,steps):sim.solve(tf,160);
}
