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
	const double s_y = 0.85e9;		// Yield stress (Pa)
	const double c0 = 0.4;
	const double tau0 = 1e-13;
	const double kappa = 4.2;
	const double Delta = 8000/TZ;
	const double Omega = 300.*22/17;		// Fix typo in 2015 JCP
	const double eps0 = 0.3;
	const double chi_inf = 900/TZ;
	const double theta = 400/TZ;
	const double Omegaeps0 = Omega*1e-30*eps0*s_y/(TZ*kB);

	// Output filenames
	const char dfn[] = "sim_sm_d", qfn[] = "sim_sm_q";

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E = 101e9;			// Young's modulus (Pa)
	const double nu = 0.35;			// Poisson ratio
	const double K = E / (3 * (1 - 2*nu)*s_y);	// Rescaled bulk modulus
	const double mu = E / (2 * (1 + nu)*s_y);	// Rescaled shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double le = 1e-2;			// Length scale (m)
	const double rho0 = 6125;			// Density (kg/m^3)
	const double sca = 2e4;
	const double t_scale = le/sqrt(mu*s_y/rho0)*sca; // Plasticity timescale (*)
	const double visc = 0.25;			// Viscous damping
	const double tmult = 0.5;			// Direct sim. timestep multiplier
	const double adapt_fac = 1e-4;		// Adaptivity factor
	const double u_bdry = 1e-7*sca;		// Wall speed (*)
    const double lamb = 1e-4;
    double tf = 1./lamb;

	// Make the output directory if it doesn't already exist
	mkdir(qs? qfn : dfn, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ, c0, tau0, kappa, Delta, Omegaeps0, chi_inf, theta, 0);

	// Initialize the simulation
	shear_sim sim(128, 64, -2, 2, -1, 1, mu, K, visc, t_scale, tmult, adapt_fac, u_bdry, lamb, &stz, qs? qfn : dfn);
	sim.init_fields(2, 600, 200);

	// Carry out the simulation using the selected simulation method
	qs? sim.solve_quasistatic(0, tf/5, 1000, 1) : sim.solve(0, tf/5, 1000);
}
