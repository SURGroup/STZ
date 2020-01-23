#ifndef BD_SIM_HH
#define BD_SIM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>

#include "vec.hh"
#include "mat.hh"
#include "fields.hh"
#include "stz_model.hh"
#include "qs_multi_bd.hh"
#include "tgmg.hh"
#include "extrap.hh"
#include "level++.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

class extra_force;

//const double bd_double_min=std::numeric_limits<double>::min();
const double bd_double_max=std::numeric_limits<double>::max();

/** \brief A class to carry out a 2D elasticity simulation. */
class bd_sim {
	public:
		/** The number of grid points in the horizontal direction. */
		const int m;
		/** The number of grid points in the vertical direction. */
		const int n;
		/** The total number of grid points. */
		const int mn;
		/** The number of staggered grid points in the horizontal direction. */
		const int me;
		/** The number of staggered grid points in the vertical direction. */
		const int ne;
		/** The total number of staggered grid points. */
		const int mne;
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
		/** The length scale for the chi diffusion. */
		double chi_len;
		/** The viscosity. */
		const double filter_stress;
		/** A time scale associated with the plastic deformation. */
		const double t_scale;
		/** A multiplier to apply to the default timestep size. */
		const double tmult;
		/** A constant used to control the adaptive timestepping of plasticity. */
		const double adapt_fac;
		/** A pointer to the class to carry out the STZ plasticty. */
		stz_dynamics *stz;
		/** The filename of the output directory. */
		const char *filename;
		/** Flags to use to determine file output. */
		const unsigned int fflags;
		/** The level set class. */
		levelset ls;
		/** A pointer to the phi values in the level set class. */
		double *phi;
		/** A pointer to the status array in the level set class. */
		int *c;
		/** A pointer to the staggered status array in the level set
		 * class. */
		int *cc;
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
		/** A boundary smoothing factor. */
		double bc_smooth;
		/** The current frame number. */
		int f_num;
		/** A debug flag to flip the spin term to the wrong sign. */
		bool wrong_spin;
		/** A flag controlling whether to apply the chi diffusion or not. */
		bool apply_chi_diffusion;
		/** A pointer to a class to apply external forces (if any). */
		extra_force *eforce;
		bd_sim(const int m_,const int n_,const double ax_,const double bx_,
			  const double ay_,const double by_,const double mu_,const double K_,
			  const double viscosity_,const double chi_len_,const double filter_stress_,
			  const double t_scale_,const double tmult_,const double adapt_fac_,
			  stz_dynamics *stz_,const char *filename_,const unsigned int fflags_);
		~bd_sim();
		void solve(double duration,int frames);
		void solve_quasistatic(double duration,int frames,int steps);
		void step_forward(double dt);
		void step_forward_quasistatic(double dt);
		void post_process(bool all_fields=true);
		void set_staggered_status();
		void direct_step(double dt);
		void advection_step(double dt);
		void projection_step(double dt);
		void init_bar(int ca);
		void write_files(int k);
		void compute_strain();
		void print_extrema(int k);
		void extrema(double *gex);
		void compute_def_rate();
		void qs_measure(double &Q,double &max_qs,bool staggered=false);
		void initialize_chi_bicubic(int mm,int nn,double *arr,double lo=0,double sca=1);
		double velocity(double x,double y);
		inline void extrapolate(bool all_fields=true) {
			ls.extrapolate_fields(e_reg);
			if(all_fields) ls.extrapolate_staggered_fields(e_st);
		}
		void set_subfield_bounds(double xmin,double xmax,double ymin,double ymax);
	private:
		/** An object containing the configuration of the linear system
		 * to be solved using the multigrid method. */
		qs_multi_bd qsm;
		/** The multigrid solver. */
		tgmg<qs_multi_bd,vec,mat> mg;
		/** Temporary storage for used during the output routine. */
		float *buf;
		ex_ref_map e_rm;
		ex_staggered e_st;
		ex_regular e_reg;
		int obnd[8];
		void set_boundaries();
		void chi_diffusion(double dt);
		double curvature(double x,double y);
		double curvature(int i,int j,double *phip);
		void output(const char *prefix,const int mode,const int sn);
		inline void rm_one_sided(double &Xd,double &Yd,double hs,c_field &f1,c_field &f2);
		inline void rm_eno2(double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3);
		inline void rmv_one_sided(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f1,c_field &f2);
		inline void rmv_eno2(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3);
		inline st_field st_one_sided(double hs,c_field &f1,c_field &f2);
		st_field st_eno2(double hs,c_field &f0,c_field &f1,c_field &f2,c_field &f3);
		inline double eno2(double p0,double p1,double p2,double p3);
		void project_stress(double &pp,double &qp,double &sp,double &taup,double phix,double phiy);
		bool bc_deriv(int i,int j,int ni,int nj,c_field &fg);
		bool bc_deriv_proj(int i,int j,int ni,int nj,c_field &fg);
		void bicubic_velocity(double x,double y,double &qu,double &qv);
		inline double phi_x(int i,int ij) {
			return (i==0?phi[ij+1]-phi[ij]:(i==m-1?phi[ij]-phi[ij-1]:(phi[ij+1]-phi[ij-1])*0.5))*xsp;
		}
		inline double phi_y(int ij) {
			return (ij<m?phi[ij+m]-phi[ij]:(ij>=mn-m?phi[ij]-phi[ij-m]:(phi[ij+m]-phi[ij-m])*0.5))*ysp;
		}
		inline c_field& diag_line(c_field *fp,int ij,bool look_more,double phic,int d1,int d2,int d2e,c_field &fg);
		void net_force(c_field *fp,int i,int j,double &fx,double &fy,double xf,double yf,c_field* fg);
		double adaptive_plastic_term(double sbar,double &chiv,double &ddev,double dt);
		inline void extrema_init(double *pex) {
			double *pe=pex+16;
			while(pex<pe) {
				*(pex++)=bd_double_max;
				*(pex++)=-bd_double_max;
			}
		}
		inline void extrema_check(double val,double *&pex) {
			if(val<*pex) *pex=val;
			pex++;
			if(val>*pex) *pex=val;
			pex++;
		}
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif
		friend class extra_force;
};

#endif
