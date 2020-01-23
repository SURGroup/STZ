#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "stz_model.hh"

const double mu=90.0/3;

double integrate(double chi0,double T,double omega,bool output) {
	const double rt3=sqrt(3.0);
	stz_dynamics_nonlinear stz(T);

	int i=0;
	const double cutoff=1e-3;
	double dx,xi=0,s=1e-20,Dp,chi=chi0,s_inf=0,dchi1,dchi2,dxm,sbar;
	bool run=true;
	while(run) {

		// Guess integration step
		dx=0.000001*(1-xi*xi);
		if(dx<1e-9) dx=1e-9;
		if(xi+dx>1-cutoff) {dx=1-cutoff-xi;run=false;}

		// Output information if requested
		sbar=abs(s)*sqrt(3);
		Dp=stz.Dplastic(sbar,chi,dchi1,dchi2);
		if(output&&++i==100) {
			printf("%g %g %g %g\n",xi,chi,s,Dp);
			i=0;
		}

		// Half of trapezoidal rule
		s_inf+=xi<1e-20?0:0.5*s*dx/xi;

		// Update fields
		xi+=0.5*dx;
		dxm=dx/(omega*xi*(1-xi*xi));
		chi+=dxm*(dchi1*(stz.chi_inf-chi)+dchi2*(stz.theta-chi));
//		chi=(chi+(dchi1*stz.chi_inf+dchi2*stz.theta)*dxm)/(1+dxm*(dchi1+dchi2));
		if(s<0) Dp=-Dp;
		s+=dx*2*mu*(xi*xi-Dp/(rt3*omega*xi))/(1-xi*xi);
		xi+=0.5*dx;

		// Half of trapezoidal rule
		s_inf+=0.5*s*dx/xi;
	}

	// Assume function is slowly varying after cutoff
	return 6*(s_inf+cutoff*s);
}

double integrate_simple(bool output) {
	int i=0;
	double dx=1/10000.0,xi=0,s_inf=0,s=1e-20;
	while(xi<1) {
		s_inf+=xi<1e-20?0:0.5*s*dx/xi;

		// Output information if requested
		if(output&&++i==10) {
			printf("%g %g\n",xi,s);
			i=0;
		}

		// Update fields
		xi+=0.5*dx;
		s+=dx*2*mu*xi/(1-xi*xi);
		//if(s>1/rt3) s=1/rt3;
		if(s>1) s=1;
		xi+=0.5*dx;

		// Half of trapezoidal rule
		s_inf+=0.5*s*dx/xi;
	}

	// Assume function is slowly varying after cutoff
	return 2*s_inf;
}

int main(int argc,char **argv) {
	if(argc!=3) {
		fputs("Syntax: ./cavity <temperature> <increase of chi0>\n",stderr);
		return 1;
	}
	double T=atof(argv[1]),chi0=(T+atof(argv[2]))/16500,omega=1e-13;

	double s_inf=integrate(chi0,T,omega/1e-13,true);
//	double s_inf=integrate_simple(true);
	fprintf(stderr,"%g\n",s_inf);
//	return 0;

//	for(double k=10;k<16;k+=0.2) {
//		omega=pow(0.1,k)/tau0;
//		printf("%g %g %g\n",tau0*omega,omega,integrate(chi0,theta_in,omega,false));
//	}
}
