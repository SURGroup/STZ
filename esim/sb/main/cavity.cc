#include <cstdio>
#include <cstdlib>
#include <cmath>

const double c0=0.4;
const double rt3=sqrt(3);
const double tau0=1e-13;
const double kappa=4.2;
const double T0=370;
const double TA=700;
const double T1=13000;
const double TZ=16500;
const double chi_inf=900/TZ;
const double Delta=0.2; // Based on 3500K
const double Omegaeps0=Delta/1.2; // Based on Omega=133 Angstrom^3 and eps0=0.3
const double alpha0=6;
const double ez=1;
const double s0=1;
const double mu=101/1.1/(2*(1+0.35));
double theta,rho_theta;

double rho(double theta) {
	const double a=2;
	double T=theta*TZ;
	return T>T0?exp(-T1/(T-T0)*exp(-a*(T-T0)/(TA-T0))-alpha0):0;
}

double sR(double sbar) {
	return exp(-Delta*exp(-Omegaeps0*sbar/Delta)/theta);
}

double sC(double sbar) {
	double arg=Omegaeps0*sbar/Delta;
	return exp(-Delta*cosh(arg)/theta)*cosh(Delta*sinh(arg)/theta);
}

double sT(double sbar) {
	return tanh(Delta*sinh(Omegaeps0*sbar/Delta)/theta);
}

double M(double sbar) {
	double x=1+sbar*sT(sbar)/s0+sC(0)*rho_theta/(2*sC(sbar));
	return s0/(2*sbar)*x-s0/(2*sbar)*sqrt(x*x-4*sbar/s0*sT(sbar));
}

double Dpl(double s,double chi) {
	double sbar=rt3*abs(s);
	return exp(-ez/chi)*sC(sbar)*(sT(sbar)-M(sbar))*s/(sbar*tau0);
}

double chich(double Dpl,double s,double chi) {
	double Gamma=6*tau0*s*Dpl/(s0*exp(-ez/chi));
	return exp(-ez/chi)*(Gamma*(chi_inf-chi)+kappa*exp(-Delta/theta)*rho_theta*(theta-chi));

}

double integrate(double chi0,double theta_in,double omega,bool output) {
	theta=theta_in;
	rho_theta=rho(theta);

	int i=0;
	const double cutoff=1e-3;
	double dx,xi=0,s=1e-20,Dp,chi=chi0,s_inf=0;
	bool run=true;
	while(run) {

		// Guess integration step, and
		dx=0.000001*(1-xi*xi);
		if(dx<1e-9) dx=1e-9;
		if(xi+dx>1-cutoff) {dx=1-cutoff-xi;run=false;}

		// Output information if requested
		Dp=Dpl(s,chi);
		if(output&&++i==100) {
			printf("%g %g %g %g\n",xi,chi,s,Dp);
			i=0;
		}

		// Half of trapezoidal rule
		s_inf+=xi<1e-20?0:0.5*s*dx/xi;

		// Update fields
		xi+=0.5*dx;
		chi+=dx*chich(Dp,s,chi)/(tau0*c0*omega*xi*(1-xi*xi));
		s+=dx*2*mu*(xi*xi-Dp/(omega*xi))/(1-xi*xi);
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
		s+=dx*2*mu*xi*xi/(1-xi*xi);
		if(s>1/rt3) s=1/rt3;
		//if(s>1) s=1;
		xi+=0.5*dx;

		// Half of trapezoidal rule
		s_inf+=0.5*s*dx/xi;
	}

	// Assume function is slowly varying after cutoff
	return 6*s_inf;
}

int main(int argc,char **argv) {
	if(argc!=3) {
		fputs("Syntax: ./cavity <temperature> <increase of chi0>\n",stderr);
		return 1;
	}
//	double theta_in=atof(argv[1])/TZ,chi0=theta_in+atof(argv[2])/TZ,omega=1e-13;

	double s_inf=integrate_simple(false);//(chi0,theta_in,omega/tau0,true);
//	double s_inf=integrate_simple(true);
	fprintf(stderr,"%g\n",s_inf);
//	return 0;

//	for(double k=10;k<16;k+=0.2) {
//		omega=pow(0.1,k)/tau0;
//		printf("%g %g %g\n",tau0*omega,omega,integrate(chi0,theta_in,omega,false));
//	}
}
