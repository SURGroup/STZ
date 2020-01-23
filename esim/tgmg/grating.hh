#ifndef TGMG_HELMHOLTZ_HH
#define TGMG_HELMHOLTZ_HH

#include <complex>
typedef std::complex<double> cpx;

#include "tgmg.hh"

struct h_grating {
	public:
		/** Grid dimensions. */
		const int m;
		const int n;
		/** Total number of gridpoints. */
		const int mn;
		/** Periodicity in the x and y directions. */
		const bool x_prd;
		const bool y_prd;
		/** The mode to use for the Gauss-Seidel smoothing. (0=default) */
		static const char gs_mode=0;
		/** Threshold on L_2 norm of residual to terminate the multigrid solve. */
		const double acc;
		/** Lower and upper limits in the x direction. */
		const double ax,bx;
		/** Lower and upper limits in the y direction. */
		const double ay,by;
		/** Grid spacings in the x and y directions. */
		const double dx,dy;
		/** A pointer to the spatially varying k field. */
		double *kf;
		/** A pointer to the memory for the linear system. */
		cpx *A;
		/** A pointer to the allocated memory. */
		cpx *b;
		/** A pointer to the solution vector. */
		cpx* x;
		cpx* z;
		h_grating(const int m_,const int n_,const double ax_,const double bx_,
			const double ay_,const double by_,const double k);
		~h_grating();
		/** Function to determine whether a grid point is on the edge or not.
		 */
		inline cpx a_dl(int i,int ij) {return A[10*ij];}
		inline cpx a_dc(int i,int ij) {return A[10*ij+1];}
		inline cpx a_dr(int i,int ij) {return A[10*ij+2];}
		inline cpx a_cl(int i,int ij) {return A[10*ij+3];}
		inline cpx a_cc(int i,int ij) {return A[10*ij+4];}
		inline cpx a_cr(int i,int ij) {return A[10*ij+5];}
		inline cpx a_ul(int i,int ij) {return A[10*ij+6];}
		inline cpx a_uc(int i,int ij) {return A[10*ij+7];}
		inline cpx a_ur(int i,int ij) {return A[10*ij+8];}
		inline cpx inv_cc(int i,int ij,cpx v) {return A[10*ij+9]*v;}
		/** Calculates the ith component of the multiplication (A-D)z, needed
		 * in the Gauss--Seidel smoothing iteration. */
		inline cpx mul_a(int i,int ij) {
			cpx *zp=z+ij,*zl=i==0?zp+(m-1):zp-1,*zr=i==m-1?zp+(1-m):zp+1,*Ap=A+10*ij;
			cpx ans;
			ans=*zl*Ap[3]+*zr*Ap[5];
			ans+=ij<m?0:Ap[0]*zl[-m]+Ap[1]*zp[-m]+Ap[2]*zr[-m];
			ans+=ij>=mn-m?0:Ap[6]*zl[m]+Ap[7]*zp[m]+Ap[8]*zr[m];
			return ans;
		}
		void solve();
		void mg_solve() {
			z=x;
			mg.verbose=2;
			mg.setup();
			mg.solve_v_cycle();
		}
		void setup_linear_system();
		void setup_test_problem();
		void output(const char* filename,cpx *ff,int mode);
	private:
		tgmg<h_grating,cpx,cpx> mg;
		cpx iprod(cpx *u,cpx *v);
		void pc_solve(cpx *u,cpx *v);
		void mul_fa(cpx *u,cpx *v);
		double l2_error(cpx *u);
};

#endif
