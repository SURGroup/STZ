#ifndef TGMG_HELMHOLTZ_HH
#define TGMG_HELMHOLTZ_HH

#include <complex>
typedef std::complex<double> cpx;
#include "tgmg.hh"

struct helmholtz {
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
		/** Stencil entries. */
		const cpx fm,fm_inv,fm_full,fex,fey;
		/** A pointer to the allocated memory. */
		cpx *b;
		/** A pointer to the solution vector. */
		cpx* x;
		cpx* z;
		helmholtz(const int m_,const int n_,const double ax_,const double bx_,
			const double ay_,const double by_,const double k);
		~helmholtz();
		/** Function to determine whether a grid point is on the edge or not.
		 */
		inline bool cent(int i,int ij) {return ij>30*m&&ij<50*m&&i>70&&i<80;}
		inline bool edge(int i,int ij) {return ij>=mn-m||ij<m;}
		/** Functions to specify the corner stencil entries. */
		inline cpx a_dl(int i,int ij) {return 0.;}
		inline cpx a_dr(int i,int ij) {return 0.;}
		inline cpx a_ul(int i,int ij) {return 0.;}
		inline cpx a_ur(int i,int ij) {return 0.;}
		/** Functions to specify the vertical stencil entries. */
		inline cpx a_dc(int i,int ij) {return edge(i,ij)?0:fey;}
		inline cpx a_uc(int i,int ij) {return edge(i,ij)?0:fey;}
		/** Functions to specify the horizontal stencil entries. */
		inline cpx a_cl(int i,int ij) {return fex;}
		inline cpx a_cr(int i,int ij) {return fex;}
		/** Function to specify the central stencil entry (on the diagonal of
		 * the linear system). */
		inline cpx a_cc(int i,int ij) {return edge(i,ij)?fm:fm;}
		/** Function to multiply by the reciprocal of the central
		 * stencil entry. This is specified as a separate function for
		 * computational efficiency. */
		inline cpx inv_cc(int i,int ij,cpx v) {return (edge(i,ij)?fm_inv:fm_inv)*v;}
		/** Calculates the ith component of the multiplication (A-D)z, needed
		 * in the Gauss--Seidel smoothing iteration. */
		inline cpx mul_a(int i,int ij) {
			return edge(i,ij)?0:fex*(z[i==m-1?ij+1-m:ij+1]+z[i==0?ij+m-1:ij-1])
				+fey*(z[ij+m]+z[ij-m]);
		}
		void solve();
		void setup_test_problem();
		void output(const char* filename,cpx *ff,int mode);
	private:
		tgmg<helmholtz,cpx,cpx> mg;
		cpx iprod(cpx *u,cpx *v);
		void pc_solve(cpx *u,cpx *v);
		void mul_fa(cpx *u,cpx *v);
		double l2_error(cpx *u);
};

#endif
