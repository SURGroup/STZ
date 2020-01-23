#ifndef LJ_SIM_HH
#define LJ_SIM_HH

#include <cstdio>

void fatal_error(const char *p,int status);
FILE* safe_fopen(const char *fn,const char *mode);

class lj_sim {
	public:
		const int mem;
		const double sx,sy;
		const double cutoffsq;
		double damp;
		double gravity;
		double vx_fix;
		double vy_fix;
		double *x;
		double *v;
		double *a;
		double *rad;
		int *ty;
		int n;
		lj_sim(int mem_,double sx_,double sy_,double cutoff_,double damp_);
		~lj_sim();
		void put(double px,double py,double pr);
		void integrate(double dt);
		void compute_acceleration();
		void fix_atoms();
		void label_atoms();
		void load_state(FILE *fp);
		inline void load_state(const char *filename) {
			FILE *fp=safe_fopen(filename,"rb");
			load_state(fp);
			fclose(fp);
		}
		void save_state(FILE *fp);
		inline void save_state(const char *filename) {
			FILE *fp=safe_fopen(filename,"wb");
			save_state(fp);
			fclose(fp);
		}
		void output_gnuplot(FILE *fp);
		inline void output_gnuplot(const char *filename) {
			FILE *fp=safe_fopen(filename,"w");
			output_gnuplot(fp);
			fclose(fp);
		}
		void output_povray(FILE *fp);
		inline void output_povray(const char *filename) {
			FILE *fp=safe_fopen(filename,"w");
			output_povray(fp);
			fclose(fp);
		}
	private:
		bool pair_force(int i,int j,double &fx,double &fy);
};

#endif
