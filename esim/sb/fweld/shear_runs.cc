#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>

#include "shear_sim.hh"

// Temperature scale that is used to non-dimensionalize temperatures
const double TZ=21000;

// The Bolztmann constant
const double kB=1.3806503e-23;

void syntax_message() {
	fputs("Syntax: ./shear_runs <case> direct    (for direct simulation)\n"
	      "        ./shear_runs <case> qs        (for quasi-static simulation)\n",stderr);
	exit(1);
}

int main(int argc,char **argv) {

	// Check command-line arguments
	if(argc!=3) syntax_message();
	bool qs;
	if(strcmp(argv[2],"direct")==0) qs=false;
	else if (strcmp(argv[2],"qs")==0) qs=true;
	else syntax_message();

	// Read in case number and check it's valid
	int run=atoi(argv[1]);
	if(run<0||run>6) {
		fputs("Case number out of range\n",stderr);
		return 1;
	}

	// STZ model parameters (based on Vitreloy 1 BMG)
	const double c0=0.4;
	const double tau0=1e-13;
	const double kappa=4.2;
	const double Delta=8000/TZ;
	const double Omega=300;
	const double eps0=0.3;
	const double chi_inf=900/TZ;
	const double theta=400/TZ;
	const double Omegaeps0=Omega*1e-30*eps0*1.1e9/(TZ*kB);

	// Output filenames
	char fn[16];
	const char* fns[]={"p","line","line2","line3","sin","low","lowc"};
	sprintf(fn,"shr_%s%s",qs?"q":"d",fns[run]);

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Parameters labeled are scaled to be larger than
	// realistic values, but allow for a direct--quasistatic comparison.
	const double scatab[7]={1e4,1e4,1e4,1e4,1,1e4,1e4};
	double sca=scatab[run];
	const double visc=0.02;			// Viscous damping
	double t_scale=4.04645058986732e-06*sca;// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	double u_bdry=1e-7*sca;			// Wall speed (*)

	// Make the output directory if it doesn't already exist
	mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the simulation
	int m=run>=5?1280:640,bx=run>=5?8:4;
	shear_sim sim(m,160,-bx,bx,-1,1,mu,K,visc,t_scale,tmult,adapt_fac,u_bdry,&stz,fn);

	// Set up the simulation fields
	int ifn[7]={0,1,1,1,2,0,0};
	const double chitab[7]={630,600,630,660,620,480,600},
	             dchitab[7]={170,200,170,140,180,320,200};
	sim.init_fields(ifn[run],chitab[run],dchitab[run]);

	// Carry out the simulation using the selected simulation method
	int frs[7]={400,400,400,400,600,500,500};
	double dur[7]={2e6,2e6,2e6,2e6,3e6,1.25e6,1.25e6};
	int nquas[7]={50,25,25,25,50,50,50};
	qs?sim.solve_quasistatic(0,dur[run]/sca,frs[run],nquas[run]):sim.solve(0,dur[run]/sca,frs[run]);

}
