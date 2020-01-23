// This is an example file for testing the multigrid code.

#include "helmholtz.hh"

const double pi=3.1415926535897932384626433832795;

int main() {
	const int m=128,n=129;
	const double ax=-8,bx=8,ay=ax,by=bx,h=(bx-ax)/m;
	const double k=2*pi/(44*h);

	helmholtz he(m,n,ax,bx,ay,by,k);
	he.setup_test_problem();
	he.solve();
	he.output("x.real",he.x,2);
	he.output("x.imag",he.x,3);
}
