#ifndef SBSIM_HH
#define SBSIM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include "common.hh"

const double pi=3.1415926535897932384626433832795;
const double tpi=2*pi;

/** \brief A class to carry out an athermal STZ simulation. */
class sbsim {
	public:
		/** The number of grid points in the horizontal direction. */
		const int m;
		/** The number of grid points in the vertical direction. */
		const int n;
		/** The total number of grid points. */
		const int mn;
		/** The lower bound in the x direction. */
		const double ax;
		/** The lower bound in the y direction. */
		const double ay;
		/** The grid spacing in the x direction. */
		const double dx;
		/** The grid spacing in the y direction. */
		const double dy;
		/** The inverse grid spacing in the x direction. */
		const double xsp;
		/** The inverse grid spacing in the y direction. */
		const double ysp;
		/** The parameter \f$\chi_\infty\f$ in the athermal STZ model.
		 */
		const double chi_inf;
		/** The parameter \f$1/c_0\f$ in the athermal STZ model. */
		const double c_0_inv;
		/** The rescaled parameter \f$\nu\f$ in the athermal STZ model. */
		const double nu;
		/** The elastic shear modulus. */
		const double mu;
		/** The elastic bulk modulus. */
		const double K;
		/** The inverse elastic shear modulus. */
		const double mu_inv;
		/** A multiplier to apply to the default timestep computation. */
		const double tmult;
		/** The viscosity. */
		const double visc;
		/** The filename of the output directory. */
		const char *filename;
		/** The horizontal velocity field. */
		double *u;
		/** The vertical velocity field. */
		double *v;
		/** The pressure field. */
		double *p;
		/** The first deviatoric stress field. */
		double *s;
		/** The second deviatoric stress field. */
		double *tau;
		/** The effective temperature field. */
		double *chi;
		/** A field in which to store the update of u. */
		double *cu;
		/** A field in which to store the update of v. */
		double *cv;
		/** A field in which to store the update of p. */
		double *cp;
		/** A field in which to store the update of s. */
		double *cs;
		/** A field in which to store the update of tau. */
		double *ctau;
		/** A field in which to store the update of chi. */
		double *cchi;
		/** The current simulation time. */
		double time;
		sbsim(int m_,int n_,double ax_,double bx_,double ay_,double by_,
		      double chi_inf_,double c_0_inv_,double nu_,double mu_,double K_,
		      double visc_,double tmult_,const char* filename);
		~sbsim();
		void solve(double t_start,double t_end,int frames);
		void step_forward(double dt);
		void init_fields();
		void print_extrema();
		void write_files(int k);
		void output(const char *suffix,double *array,const int sn);
	private:
		/** Temporary storage used during file output. */
		float *buf;
		/** A pointer to the end of the output buffer. */
		float *bufe;
		inline double qs(double s);
		inline double eno2(double p0,double p1,double p2,double p3);
		void eno2(int ij0,int ij1,int ij2,int ij3,double &uv,double &vv,double &pv,double &sv,double &tauv,double &chiv);
};

#endif
