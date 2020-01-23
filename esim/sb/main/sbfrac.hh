#ifndef SBFRAC_HH
#define SBFRAC_HH

#include <cmath>
#include <limits>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "common.hh"
#include "fileinfo.hh"
#include "bd_sim.hh"
#include "stz_model.hh"

class sbfrac : public fileinfo {
	public:
		int me;
		bd_sim *fs;
		FILE *kfile,*ffile,*pfile;
		bool nucleated;
		sbfrac(const char *fn) : fileinfo(fn), me(m+1), nucleated(false), d_counter(0) {}
		~sbfrac() {}
		void solve();
		void init_fields();
		void print_extrema();
		void nucleate_void();
		void def_rate_computation();
		void strain_computation();
		void step_forward(double dt);
		void step_forward_quasistatic(double dt);
		void set_boundaries();
		double adaptive_ts_select();
		void write_files(int k);
		void post_process(bool all_fields=true);
	private:
		int d_counter;
		void random_chi_init(double *&rnoise,int &k,double &fac,double &scale);
		double arg(double x,double y);
		void set_irwin_velocity(c_field &f,double xs,double ys);

		inline double rshift() {
			return void_nucl_rand*(0.5-double(rand())/RAND_MAX);
		}
		inline int find_gaussian_cutoff(double fac);
		inline void box_muller(double &a,double &b) {
			double theta=2*3.1415926535897932384626433832795*rnd(),r=sqrt(-2*log(rnd()+0.5/RAND_MAX));
			a=r*cos(theta);
			b=r*sin(theta);
		}
		inline double rnd() {
			return (double(rand())+0.5)/RAND_MAX;
		}
#ifdef _OPENMP
		inline double wtime() {return omp_get_wtime();}
#else
		inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif
};

#endif
