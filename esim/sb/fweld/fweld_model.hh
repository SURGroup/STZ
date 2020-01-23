#ifndef FWELD_MODEL_HH
#define FWELD_MODEL_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

class fweld_dynamics {
	public:
		const double d0;
		fweld_dynamics(double d0_) : d0(d0_) {}
		~fweld_dynamics() {}
		double Dplastic(double sbar,double T,double &dT) {
			dT=0.01;
			return d0*sbar;
		}
};

#endif
