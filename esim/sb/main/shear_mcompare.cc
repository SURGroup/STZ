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

int main() {
	int i,j,k;

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

	// Elasticity related parameters (based on Vitreloy 1 BMG)
	const double E=101;			// Young's modulus
	const double nu=0.35;			// Poisson ratio
	const double s_y=0.85;			// Yield stress
	const double K=E/(3*(1-2*nu)*s_y);	// Bulk modulus
	const double mu=E/(2*(1+nu)*s_y);	// Shear modulus

	// Other parameters. Note that parameters labeled (*) are far larger
	// than realistic values, but allow for a direct--quasistatic
	// comparison.
	const double sca=1.25e3;
	const double visc=0.02;			// Viscous damping
	const double t_scale=4.04645058986732e-06*sca;// Plasticity timescale (*)
	const double tmult=0.5;			// Direct sim. timestep multiplier
	const double adapt_fac=0.002;		// Adaptivity factor
	const double u_bdry=1e-7*sca;		// Wall speed (*)

	// Initialize an STZ plasticity model
	stz_dynamics_linear_athermal stz(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,0);

	// Initialize the direct simulation
	shear_sim dsim(640,160,-4,4,-1,1,mu,K,visc,0,0,t_scale,tmult,adapt_fac,u_bdry,&stz,"");
	dsim.init_fields(1,600,200);
	dsim.time=0;

	// Initialize the quasi-static simulations
	shear_sim* qsim[7];
	for(i=0;i<7;i++) {
		qsim[i]=new shear_sim(640,160,-4,4,-1,1,mu,K,visc,0,0,t_scale,tmult,adapt_fac,u_bdry,&stz,"");
		qsim[i]->init_fields(1,600,200);
		qsim[i]->time=0;
	}

	// Set up constants for the simulations
	int ti[7]={5,10,20,40,80,160,320};
	double l2[4],target_time,time_interval=1000/sca,qdt[7];
	const double dt=dsim.dx*dsim.dx*dsim.tmult;

	// Open output file
	FILE *fp=fopen("m_comp","w");
	if(fp==NULL) {
		fputs("Can't open output file\n",stderr);
		return 1;
	}

	// Carry out first comparison
	fprintf(fp,"%g",dsim.time);
	for(i=0;i<7;i++) {
		qdt[i]=time_interval/ti[i];
		dsim.l2_comparison(*qsim[i],l2);
		fprintf(fp," %.12g %.12g %.12g %.12g",*l2,l2[1],l2[2],l2[3]);
	}
	fputs("\n",fp);

	// Advance the simulations forward concurrently
	for(i=1;i<=2000;i++) {
		fflush(fp);

		// Direct simulation update
		target_time=time_interval*i;
		while(dsim.time+dt*(1+1e-8)<target_time) dsim.step_forward(dt);
		dsim.step_forward(target_time-dsim.time);

		// Quasi-static simulation update and output
		fprintf(fp,"%g",dsim.time);
		for(j=0;j<7;j++) {
			for(k=0;k<ti[j];k++) {printf("%d %d ",j,k);qsim[j]->step_forward_quasistatic(qdt[j]);}
			dsim.l2_comparison(*qsim[j],l2);
			fprintf(fp," %.12g %.12g %.12g %.12g",*l2,l2[1],l2[2],l2[3]);
		}
		fputs("\n",fp);
	}
	fclose(fp);

	// Free memory for dynamically allocated shear classes
	for(i=0;i<7;i++) delete qsim[i];
}
