#include <cstdio>

// Lowest interpolation data point
const double lo=150.;

// Step between interpolation data points, and inverse step
const double step=100.;
const double istep=1./step;

// Number of data points, and data values
const int n=8;
const double val[n]={13.,12.5,11.2,8.5,5.5,3,2,1.8};

// Computes the linear interpolation for a given function value
double lin_interp(double x) {
	x=(x-lo)*istep;
	int i=int(x);
	if(i<0) i=0;else if(i>n-2) i=n-2;
	x-=i;
	return val[i]*(1-x)+val[i+1]*x;
}

int main() {

	// Loop over a sequence of x values and compute the corresponding function value
	for(double x=0.;x<1001.;x+=10.)
		printf("%g %g\n",x,lin_interp(x));
}
