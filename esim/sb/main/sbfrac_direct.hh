#ifndef SBFRAC_HH
#define SBFRAC_HH

#include <cstdlib>

const double tpi=2*pi;

class sbfrac {
	public:
		const int m,n,mn;
		const double ax,ay,dx,dy,xsp,ysp;
		const double chi_inf,c_0_inv,nu,mu,K,mu_inv,tmult;
		const char *filename;
		ofstream kfile;
		levelset ls;
		double *phi;int *c;
		double *u,*v,*p,*s,*tau,*chi;
		double *cu,*cv,*cp,*cs,*ctau,*cchi;
#ifdef TRACERS
		double *trace;
#endif
		double time;
		sbfrac(int im,int in,double iax,double ibx,double iay,double iby,
		      double ichi_inf,double ic_0_inv,double inu,double imu,double iK,double itimestep,
		      const char* filename);
		~sbfrac();
		void solve(double t_start,double t_end,int frames);
		void init_fields();
		void print_maxima();
		void write_files(int k);
		bool compute_bc_deriv(double &uu,double &vu,double &pu,double &su,double &tauu,double &chiu,int i,int j,int ni,int nj);
		template <class T>
		void output(const char *filename,const char *suffix,T *array,const int sn);
		void clean();
		void nucleate_void();
		void relax_fields();
		double velocity(double x,double y);
	private:
		void fatal_error(const char *p,int code);
		inline double phi_x(int i,int ij);
		inline double phi_y(int i,int ij);
		inline double qs(double s);
		inline double qs_new(double s,double cchi);
		inline double eno2(double p0,double p1,double p2,double p3);
		inline double arg(double x,double y);
		inline double frho(double T);
		inline double fsC(double s);
};

#endif
