#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>

using std::string;

#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures.
const double TZ = 21000;

// The Boltzmann constant.
const double kB = 1.3806503e-23;

int main(int argc, char **argv) {
    // STZ model parameters (based on Vitreloy 1 BMG).
    const double c0        = 0.4;
    const double tau0      = 1e-13;
    const double kappa     = 4.2;
    const double Delta     = 8000/TZ;
    const double Omega     = 300.*22./17.;
    const double eps0      = 0.3;
    const double chi_inf   = 900/TZ;
    const double theta     = 400/TZ;
    const double Omegaeps0 = Omega*1e-30*eps0*1.1e9/(TZ*kB);

    // Elasticity related parameters (based on Vitreloy 1 BMG).
    const double E   = 101;                     // Young's modulus
    const double nu  = 0.35;                    // Poisson ratio
    const double s_y = 0.85;                    // Yield stress
    const double K   = E/(3*(1 - 2*nu)*s_y);    // Bulk modulus
    const double mu  = E/(2*(1 + nu)*s_y);      // Shear modulus

    // Scaling factor from physical values.
    double sca = 1e4;

    // Other parameters.
    const double visc      = 0.25;              // Viscous damping
    const double chi_len   = 0.0;               // Dpl-mediated diffusion
    const double le        = 1e-2;
    const double rho0      = 6125;
    double t_scale         = le/sqrt(mu*s_y/rho0)*sca; // Plasticity timescale (*)
    const double tmult     = 0.5;               // Direct sim. timestep multiplier
    const double adapt_fac = 1e-4;              // Adaptivity factor
    const double lamb      = 1e-2;
    const double u_bdry    = lamb;

    // Variables for advancing the two simulations forward in time,
    // and storing the l2 comparison result.
    double l2[4];

    // Initialize an STZ plasticity model.
    stz_dynamics_linear_athermal stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);

    // Initialize the two simulations.
    // ss contains different numbers of grid points for l2 comparison between transformed and untransformed simulations.
    const int ss[9] = {64, 128, 256, 384, 512, 640, 768, 896, 1024};

    // Output fields and miscellaneous flags. 1-u, 2-v, 4-p, 8-q, 16-s, 32-tau,
    // 64-chi, 128-tem, 256-dev, 512-X, 1024-Y, 2048-(qs measure),
    // 4096-(wtf check).
    const unsigned int fflags=4096;

    // Timestep of the first simulation, and then define a final time based on that value.
    double tf, dt;

    // Open the output file.
     char buf[16];
    sprintf(buf, "sc_test.dat");
    FILE *fp = fopen(buf, "w");
    if (fp == NULL) {
        fputs("Can't open output file\n", stderr);
        return 1;
    }

    // Loop over all the discretizations.
    for (int curr_run = 0; curr_run < 9; curr_run++) {
        // Step both simulations forward until they hit tf, so that we can compare
        // the simulation data at a fixed time point.
        fflush(fp);

        // Dynamically allocate memory for the current simulations.
        shear_sim *t_ptr = new shear_sim(2*ss[curr_run], ss[curr_run], -2, 2, -1, 1, mu, K, visc, chi_len, t_scale, adapt_fac,      0, lamb, &stz, false, fflags, "");

        shear_sim *p_ptr = new shear_sim(2*ss[curr_run], ss[curr_run], -2, 2, -1, 1, mu, K, visc, chi_len, t_scale, adapt_fac, u_bdry,    0, &stz, false, fflags, "");

        // Dereference the pointers for handles on the current simulations.
        shear_sim &trans_sim = *t_ptr;
        shear_sim &phys_sim  = *p_ptr;

        // Initialize the field values.
        trans_sim.time = phys_sim.time = 0;
        trans_sim.init_fields(1, 600, 200);
        phys_sim.init_fields(1, 600, 200);

        // Update the current value of dt.
        dt = trans_sim.dx*trans_sim.dx*tmult;

        // Define the final time.
        if (curr_run == 0) {
            tf = 1;
        }

        // Transformed simulation update.
        while (trans_sim.time + dt*(1 + 1e-11) < tf) trans_sim.step_forward(dt);
        trans_sim.step_forward(tf - trans_sim.time);

        // Physical simulation update.
        while (phys_sim.time + dt*(1 + 1e-11) < tf) phys_sim.step_forward(dt);
        phys_sim.step_forward(tf - phys_sim.time);

        // Compare fields.
        trans_sim.l2_comparison_transform(phys_sim, l2);
        fprintf(fp, "%g %.12g %.12g %.12g %.12g\n", trans_sim.time, *l2, l2[1], l2[2], l2[3]);

        delete t_ptr;
        delete p_ptr;
    }

    // Close the output file.
    fclose(fp);
}
