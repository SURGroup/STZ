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
    fputs("Syntax: ./shear_small direct    (for direct simulation)\n"
          "        ./shear_small qs        (for quasi-static simulation)\n",stderr);
    exit(1);
}

int main(int argc,char **argv) {

    // Check command-line arguments
    if (argc != 2) syntax_message();
    bool qs = true;
    if (strcmp(argv[1], "direct") == 0) qs = false;
    else if (strcmp(argv[1], "qs") != 0) syntax_message();

    // STZ model parameters (based on Vitreloy 1 BMG)
    const double s_y = 0.85e9;              // Yield stress (Pa)
    const double c0 = 0.4;
    const double tau0 = 1e-13;
    const double kappa = 4.2;
    const double Delta = 8000/TZ;
    const double Omega = 300.*22/17;        // Fix typo in 2015 JCP
    const double eps0 = 0.3;
    const double chi_inf = 900/TZ;
    const double theta = 400/TZ;
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
    const double le = 1e-2;                 // Length scale (m)
    const double rho0 = 6125;               // Density (kg/m^3)
    const double sca = 2e4;
    const double t_scale = le/sqrt(mu*s_y/rho0)*sca; // Plasticity timescale (*)
    const double visc = 0.02;               // Viscous damping
    const double chi_len=0.1;               // Dpl-mediated diffusion
    const double adapt_fac = 1e-4;          // Adaptivity factor

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
    stz_dynamics_linear_athermal stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);

    // Initialize the simulation
    shear_sim sim(200, 100, -2, 2, -1, 1, mu, K, visc, chi_len, t_scale, adapt_fac, u_bdry, lamb, &stz, y_prd, fflags, qs? qfn : dfn);
    sim.init_fields(0, 580, 220);
    sim.initialize_tracers(40,20);

    // Set up an example input field to use for bicubic interpolation
    // of the effective temperature
    const int mm=64,nn=32;
    double f[mm*nn];
    for(int j=0;j<nn;j++) {
        double *fp=f+j*mm;
        for(int i=0;i<mm;i++,fp++) *fp=i>30&&i<33&&j>14&&j<17;
    }
    sim.initialize_chi_bicubic(mm,nn,f,550./TZ,30./TZ);

    // Carry out the simulation using the selected simulation method
    qs?sim.solve_quasistatic(tf,120,5):sim.solve(tf,160);
}
