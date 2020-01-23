#ifndef SBSIM_HH
#define SBSIM_HH

const double pi=3.1415926535897932384626433832795;
const double tpi=2*pi;

class sbsim {
	public:
		const int m,n,mn;
		const double ax,ay,dx,dy,xsp,ysp;
		const double chi_inf,c_0_inv,nu,mu,K,mu_inv,tmult;
		const char *filename;
		std::ofstream wallfile;
		levelset ls;
		double *phi;int *c;
		double *u,*v,*p,*s,*tau,*chi;
		double *cu,*cv,*cp,*cs,*ctau,*cchi;
		int min_i,max_i;
		double wallx,wallu,wallv,time;
		sbsim(int im,int in,double iax,double ibx,double iay,double iby,
		      double ichi_inf,double ic_0_inv,double inu,double imu,double iK,double itimestep,
		      const char* filename);
		~sbsim();
		void solve(double t_start,double t_end,int frames);
		void init_fields();
		void print_maxima();
		void write_files(int k);
		bool compute_bc_deriv(double &uu,double &vu,double &pu,double &su,double &tauu,double &chiu,int i,int j,int ni,int nj);
		template <class T>
		void output(const char *filename,const char *suffix,T *array,const int sn);
		void clean();
		void nucleate_void();
		double velocity(double x,double y);
	private:
		ex_array e_reg;
		void fatal_error(const char *p,int code);
		inline double phi_x(int i,int ij);
		inline double phi_y(int i,int ij);
		inline double qs(double s);
		inline double qs_new(double s,double cchi);
		inline double eno2(double p0,double p1,double p2,double p3);
		void simple_extrapolation(int i,int i2,int i3);
};

#endif
