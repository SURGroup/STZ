#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "sbfrac.hh"

int main(int argc,char **argv) {
	if(argc!=2) {
		fprintf(stderr,"Usage: ./crack_tip <sim_name>\n");
		return 1;
	}

	sbfrac sim(argv[1]);
	sim.solve();
}
