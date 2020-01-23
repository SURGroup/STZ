// This is an example file for testing the multigrid code.

#include "grating.hh"

const double pi=3.1415926535897932384626433832795;

int main() {
	const int m=128,n=129;
	const double ax=-0.5,bx=0.5,ay=0,by=1;
	const double k=1;

	h_grating he(m,n,ax,bx,ay,by,k);
	he.setup_test_problem();
	he.setup_linear_system();
	he.mg_solve();
	he.output("x.real",he.x,2);
	he.output("x.imag",he.x,3);
}
