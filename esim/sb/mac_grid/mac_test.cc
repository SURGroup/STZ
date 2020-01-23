#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "mac_sim.hh"

const char fn[]="mtest.out";

int main(int argc,char **argv) {

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=11;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=1.0;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters.
	const double visc=0.005;		// Viscous damping
	const double tmult=0.5;			// Direct sim. timestep multiplier

	// Make the output directory if it doesn't already exist
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize the simulation
	mac_sim sim(200,200,-2,2,-2,2,mu,K,visc,tmult,fn);
	sim.init_fields();

	// Carry out the simulation using the selected simulation method
	sim.solve(4,400);
}
