#include <cstdio>

#include "fweld_model.hh"

int main() {
	fweld_dynamics fwd(2.0);
	const double T=1000.;
	double dT,val;

	for(double sbar=0;sbar<4.0;sbar+=0.1) {
		val=fwd.Dplastic(sbar,T,dT);
		printf("%g %g %g\n",sbar,val,dT);
	}
}
