#ifndef MAC_SIM_HH
#define MAC_SIM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
#include "omp.h"
#else
#include <ctime>
#endif

#include "fields.hh"

/** \brief A class to carry out a 2D elasticity simulation. */
class mac_sim {
	public:
		/** The number of grid points in the horizontal direction. */
		const int m;
		/** The number of grid points in the vertical direction. */
		const int n;
		/** The total number of grid points. */
		const int mn;
		/** The memory step length, taking into account ghost point allocation. */
		const int ml;
		/** The lower bound in the x direction. */
		const double ax;
		/** The lower bound in the y direction. */
		const double ay;
		/** The upper bound in the x direction. */
		const double bx;
		/** The upper bound in the y direction. */
		const double by;
		/** The grid spacing in the x direction. */
		const double dx;
		/** The grid spacing in the y direction. */
		const double dy;
		/** The inverse grid spacing in the x direction. */
		const double xsp;
		/** The inverse grid spacing in the y direction. */
		const double ysp;
		/** The elastic shear modulus. */
		const double mu;
		/** The reciprocal of the elastic shear modulus. */
		const double mu_inv;
		/** The elastic bulk modulus. */
		const double K;
		/** The viscosity. */
		const double viscosity;
		/** A multiplier to apply to the default timestep size. */
		const double tmult;
		/** The filename of the output directory. */
		const char *filename;
		/** An array containing the simulation fields. */
		field* const fbase;
		/** A pointer to the (0,0) grid cell in the field array. */
		field* const fm;
		/** The current simulation time. */
		double time;
		/** The current frame number. */
		int f_num;
		mac_sim(const int m_,const int n_,const double ax_,const double bx_,
			  const double ay_,const double by_,const double mu_,const double K_,
			  const double viscosity_,const double tmult_,const char *filename_);
		~mac_sim();
		void solve(double duration,int frames);
		void step_forward(double dt);
		void init_fields();
		void write_files(int k);
	private:
		void set_boundaries();
		void output(const char *prefix,const int mode,const int sn,const bool ghost=true);
		inline void net_force(field *fp,double xf,double yf,double &fx,double &fy);
		inline void c_eno2(double &uv,double &vv,double &s33v,double &Xv,double &Yv,double hs,field &f0,field &f1,field &f2,field &f3);
		inline void l_eno2(double &s11v,double &s21v,double hs,field &f0,field &f1,field &f2,field &f3);
		inline void d_eno2(double &s21v,double &s22v,double hs,field &f0,field &f1,field &f2,field &f3);
		inline double eno2(double p0,double p1,double p2,double p3);
		double pressure(field *fp);
		double dev_sq(field *fp);
		inline double dev(field *fp) {
			return sqrt(dev_sq(fp));
		}
		/** Temporary storage for used during the output routine. */
		float *buf;
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return double(clock())/CLOCKS_PER_SEC;}
#endif
};

#endif
