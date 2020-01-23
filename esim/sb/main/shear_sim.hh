#ifndef SHEAR_SIM_HH
#define SHEAR_SIM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "vec.hh"
#include "mat.hh"
#include "fields.hh"
#include "stz_model.hh"
#include "qs_multi.hh"
#include "tgmg.hh"

/** The number of standard deviations to consider when doing Gaussian
 * smoothing. */
const double gaussian_fac=5.;

/** \brief A class to carry out a 2D elasticity simulation. */
class shear_sim {
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
		/** The upper bound in the x direction. */
		const double bx;
		/** The upper bound in the y direction. */
		const double by;
		/** The simulation width. */
		const double lx;
		/** The simulation height. */
		const double ly;
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
		/** The amount of regular effective temperature diffusion. */
		const double chi_diff;
		/** The amount of Dpl-mediated effective temperature diffusion.
		 */
		const double chi_len;
		/** A measure of what one unit of simulation time corresponds
		 * to in seconds. This is used to properly scale the plastic
		 * deformation. */
		const double t_scale;
		/** A multiplier to apply to the default timestep size. */
		const double tmult;
		/** A constant used to control the adaptive timestepping of
		 * plasticity. */
		const double adapt_fac;
		/** The velocity of the top boundary. */
		const double u_bdry;
		/** The boundary condition multiplier. */
		const double bcs;
		/** A pointer to the class to carry out the STZ plasticty. */
		stz_dynamics *stz;
		/** The filename of the output directory. */
		const char *filename;
		/** An array containing the simulation fields. */
		c_field* const fm;
		/** An array for updating the velocity during the algebraic
		 * multigrid solve. */
		vec* const vel;
		/** An array for the source terms used during the algebraic
		 * multigrid solve. */
		vec* const src;
		/** The current simulation time. */
		double time;
		/** The value of T_Z, used to convert the dimensionless effective
		 * temperatures into Kelvins. */
		double TZ;
		/** The limiting value of effective temperature. */
		double chi_inf;
		/** The current frame number. */
		int f_num;
		shear_sim(const int m_,const int n_,const double ax_,const double bx_,
			  const double ay_,const double by_,const double mu_,const double K_,
			  const double visc_,const double chi_diff_,const double chi_len_,const double t_scale_,
			  const double tmult_,const double adapt_fac,const double u_bdry_,
			  stz_dynamics *stz_,const char *filename_);
		~shear_sim();
		void solve(double duration,int frames);
		void solve_quasistatic(double duration,int frames,int steps);
		void step_forward(double dt);
		void step_forward_quasistatic(double dt);
		void chi_diffusion(double dt);
		void advection_step(double dt);
		void projection_step(double dt);
		void init_fields(int chi_case,double tem_base=600,double tem_delta=200);
		void write_files(int k);
		void compute_strain();
		void qs_measure(double &Q,double &max_qs,bool staggered=false);
		void tractions(double &bx,double &by,double &tx,double &ty);
		void l2_comparison(shear_sim &ss,double *l2);
		void initialize_chi_bicubic(int mm,int nn,double *arr,double lo=0,double sca=1);
		void initialize_random(double l,double tem_avg,double tem_stdev);
		inline double avg_iters() {return qsm.tp.avg_iters();}
		/** Returns the velocity of the boundary, based on the current
		 * simulation time.
		 * \return The velocity. */
		inline double bdry_vel() {return time<1?time*u_bdry:u_bdry;}
		/** Returns the displacement of the boundary, based on the
		 * current simulation time.
		 * \return The displacement. */
		inline double bdry_pos() {return time<1?u_bdry*time*time*0.5:u_bdry*(time-0.5);}
	private:
		/** An object containing the configuration of the linear system
		 * to be solved using the multigrid method. */
		qs_multi qsm;
		void set_boundaries();
		void output(const char *prefix,const int mode,const int sn);
		inline void net_force(double &fx,double &fy,double xf,double yf,c_field *fp2,c_field *fp3);
		inline void rm_one_sided(double &Xd,double &Yd,double hs,c_field &f1,c_field &f2);
		inline void rm_eno2(double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3,double X0=0,double X2=0,double X3=0);
		inline void rmv_one_sided(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f1,c_field &f2);
		inline void rmv_eno2(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3,double X0=0,double X2=0,double X3=0);
		inline st_field st_one_sided(double hs,c_field &f1,c_field &f2);
		inline st_field st_eno2(double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3);
		inline double eno2(double p0,double p1,double p2,double p3);
		inline void net_force(double &fx,double &fy,double xf,double yf,c_field *fp,int d);
		double adaptive_plastic_term(double sbar,double &chiv,double &ddev,double dt);
		double gaussian_normal_factor(double llinv,int cut);
		void box_muller(double &r0,double &r1);
		inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
		/** Temporary storage for used during the output routine. */
		float *buf;
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return 0;}
#endif
};

#endif
