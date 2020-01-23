// This is an example file for testing the multigrid code.

#include "tgmg.cc"

// Multisetup structure for a Poisson problem with periodicity
struct multisetup1 {
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
	/** Lower and upper limits in the x direction. */
	const double ax,bx;
	/** Lower and upper limits in the y direction. */
	const double ay,by;
	/** Grid spacings in the x and y directions. */
	const double dx,dy;
	/** Ratios of grid spacings (used in some finite-element formulations). */
	const double dxdy,dydx;
	/** Stencil entries. */
	const double fm,fm_inv,fex,fey;
	/** Threshold on L_2 norm of residual to terminate the multigrid solve. */
	const double acc;
	/** A pointer to the solution vector. */
	double* const z;
	multisetup1(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_,double* const z_)
		: m(m_), n(n_), mn(m_*n_), x_prd(true), y_prd(false), ax(ax_),
		bx(bx_), ay(ay_), by(by_), dx((bx-ax)/(m-1)),
		dy((by-ay)/(n-1)), dxdy(dx/dy), dydx(dy/dx), fm(-4/(dx*dx)),
		fm_inv(-0.25*dx*dx), fex(1/(dx*dx)), fey(1/(dx*dx)),
		acc(tgmg_accuracy(fm,1e4)), z(z_) {}
	/** Function to determine whether a grid point is on the edge or not.
	 */
	inline bool edge(int i,int ij) {return (!x_prd&&(i==0||i==m-1))||(!y_prd&&(ij>=mn-m||ij<m));}
	/** Functions to specify the corner stencil entries. */
	inline double a_dl(int i,int ij) {return 0;}
	inline double a_dr(int i,int ij) {return 0;}
	inline double a_ul(int i,int ij) {return 0;}
	inline double a_ur(int i,int ij) {return 0;}
	/** Functions to specify the vertical stencil entries. */
	inline double a_dc(int i,int ij) {return edge(i,ij)?0:fey;}
	inline double a_uc(int i,int ij) {return edge(i,ij)?0:fey;}
	/** Functions to specify the horizontal stencil entries. */
	inline double a_cl(int i,int ij) {return edge(i,ij)?0:fex;}
	inline double a_cr(int i,int ij) {return edge(i,ij)?0:fex;}
	/** Function to specify the central stencil entry (on the diagonal of
	 * the linear system). */
	inline double a_cc(int i,int ij) {return fm;}
	/** Function to multiply by the reciprocal of the central stencil
	 * entry. This is specified as a separate function for computational
	 * efficiency. */
	inline double inv_cc(int i,int ij,double v) {return fm_inv*v;}
	/** Calculates the ith component of the multiplication (A-D)z, needed
	 * in the Gauss--Seidel smoothing iteration. */
	inline double mul_a(int i,int ij) {
		if(edge(i,ij)) return 0;
		int l=-1,r=1,u=m,d=-m;
		if(x_prd) {if(i==m-1) r-=m;else if(i==0) l+=m;}
		if(y_prd) {if(ij>=mn-m) u-=mn;else if(ij<m) d+=mn;}
		double *zp=z+ij;
		return fex*(zp[l]+zp[r])+fey*(zp[d]+zp[u]);
	}
};

int main() {
	const int m=64,n=65,mn=m*n;
	const double ax=-1+1.0/m,bx=1-1.0/m,ay=-1,by=1;
	int i,j,ij;
	double *b=new double[mn],*z=new double[mn];
	double x,y;
	multisetup1 msu(m,n,ax,bx,ay,by,z);
	tgmg<multisetup1,double,double> mg(msu,b,z);
	mg.verbose=2;

	// Set up solution and source arrays
	for(ij=j=0,y=ay;j<n;j++,y+=msu.dy) {
		for(i=0,x=ax;i<m;i++,x+=msu.dx,ij++) {
			z[ij]=0;
			b[ij]=0;
		}
	}
	b[32*64+5]=1;

	// Set up the multigrid hierarchy
	mg.setup();

	// Solve using multigrid V-cycles
	mg.solve_v_cycle();

	// Output the solutions in a format that can be read by Gnuplot using
	// the command "splot 'filename' matrix binary"
	mg.output_b("b.0");
	mg.output_z("z.0");
	mg.output_res("r.0");

	// Delete dynamically allocated memory
	delete [] z;
	delete [] b;
}
