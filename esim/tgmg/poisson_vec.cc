// This is an example file for testing the multigrid code.

#include "vec.hh"
#include "tgmg.cc"

// Multisetup structure for Poisson problem
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
	/** Threshold on L_2 norm of residual to terminate the multigrid solve. */
	const double acc;
	/** Lower and upper limits in the x direction. */
	const double ax,bx;
	/** Lower and upper limits in the y direction. */
	const double ay,by;
	/** Grid spacings in the x and y directions. */
	const double dx,dy;
	/** Ratios of grid spacings (used in some finite-element formulations). */
	const double dxdy,dydx;
	/** Stencil entries. */
	const double fm,fm_inv,fex,fey,fc;
	/** A pointer to the solution vector. */
	vec* const z;
	multisetup1(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_,vec* const z_)
		: m(m_), n(n_), mn(m_*n_), x_prd(false), y_prd(false),
		acc(1e-20), ax(ax_), bx(bx_), ay(ay_), by(by_),
		dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), dxdy(dx/dy), dydx(dy/dx),
		fm(-4/(dx*dx)), fm_inv(-0.25*dx*dx), fex(1/(dx*dx)),
		fey(1/(dx*dx)), fc(0), z(z_) {}
	/** Function to determine whether a grid point is on the edge or not.
	 */
	inline bool edge(int i,int ij) {return i==0||i==m-1||ij>=mn-m||ij<m;}
	/** Functions to specify the corner stencil entries. */
	inline double a_dl(int i,int ij) {return edge(i,ij)?0:fc;}
	inline double a_dr(int i,int ij) {return edge(i,ij)?0:fc;}
	inline double a_ul(int i,int ij) {return edge(i,ij)?0:fc;}
	inline double a_ur(int i,int ij) {return edge(i,ij)?0:fc;}
	/** Functions to specify the vertical stencil entries. */
	inline double a_dc(int i,int ij) {return edge(i,ij)?0:fey;}
	inline double a_uc(int i,int ij) {return edge(i,ij)?0:fey;}
	/** Functions to specify the horizontal stencil entries. */
	inline double a_cl(int i,int ij) {return edge(i,ij)?0:fex;}
	inline double a_cr(int i,int ij) {return edge(i,ij)?0:fex;}
	/** Function to specify the central stencil entry (on the diagonal of
	 * the linear system). */
	inline double a_cc(int i,int ij) {return edge(i,ij)?fm:fm;}
	/** Function to multiply by the reciprocal of the central stencil
	 * entry. This is specified as a separate function for computational
	 * efficiency. */
	inline vec inv_cc(int i,int ij,vec v) {return (edge(i,ij)?fm_inv:fm_inv)*v;}
	/** Calculates the ith component of the multiplication (A-D)z, needed
	 * in the Gauss--Seidel smoothing iteration. */
	inline vec mul_a(int i,int ij) {
		return edge(i,ij)?vec(0,0):fex*(z[ij+1]+z[ij-1])+fey*(z[ij+m]+z[ij-m]);
	}
};

int main() {
	const int m=1025,n=1025,mn=m*n;
	const double ax=-8,bx=8,ay=ax,by=bx;
	int i,j,ij;
	vec *b=new vec[mn],*z=new vec[mn];
	double x,y;
	multisetup1 msu(m,n,ax,bx,ay,by,z);
	tgmg<multisetup1,vec,double> mg(msu,b,z);
	mg.verbose=1;

	// Set up the multigrid hierarchy
	mg.setup();

	for(int k=0;k<10;k++) {

		// Set up the solution and source arrays
		for(ij=j=0,y=ay;j<n;j++,y+=msu.dy) {
			for(i=0,x=ax;i<m;i++,x+=msu.dx,ij++) {
				z[ij]=vec(0,0);
				b[ij]=msu.edge(i,ij)?vec(0,0):vec(1./(x*x+y*y+4),1./(x*x+y*y+8));
			}
		}

		// Solve using multigrid V-cycles
		mg.solve_v_cycle();
	}

	// Output the solutions in a format that can be read by Gnuplot using
	// the command "splot 'filename' matrix binary"
	//mg.output_b("b.0");
	//mg.output_z("z.0");
	//mg.output_res("r.0");

	// Delete dynamically allocated memory
	delete [] z;
	delete [] b;
}
