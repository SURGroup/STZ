#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

#include "common.hh"
#include "sbsim.hh"

int main() {
	const char fname[]="gauss_push";

	// Make the output directory if it doesn't already exist
	mkdir(fname,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize the simulation
	sbsim sim(150,150,-2,2,-2,2,0.15,1,1e9,30,60,0.002,0.5,"gauss_push");

	// Carry out the simulation
	sim.solve(0,10,250);
}

/** Initializes the simulation fields. */
void sbsim::init_fields() {
	int i,j,ij;
	double x,y;
	for(ij=j=0;j<n;j++) {
		y=ay+(j+0.5)*dy;
		for(i=0;i<m;i++,ij++) {
			x=ax+(i+0.5)*dx;
//			u[ij]=0;
			u[ij]=0.1*exp(-10*((x-1)*(x-1)+y*y));
			v[ij]=0;
			p[ij]=0;
			s[ij]=0;
			tau[ij]=0;
			chi[ij]=0.074;
		}
	}
}
