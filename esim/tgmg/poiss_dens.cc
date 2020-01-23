// This is an example file for testing the multigrid code.

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
	double* const z;
	/** A pointer to the FD stencil info. */
	double *e;
	multisetup1(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_,double* const z_)
		: m(m_), n(n_), mn(m_*n_), x_prd(false), y_prd(false),
		acc(1e-20), ax(ax_), bx(bx_), ay(ay_), by(by_),
		dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), dxdy(dx/dy), dydx(dy/dx),
		fm(-4/(dx*dx)), fm_inv(-0.25*dx*dx), fex(1/(dx*dx)),
		fey(1/(dy*dy)), fc(0), z(z_), e(new double[6*mn]) {}
	~multisetup1() {
		delete [] e;
	}
	void create_table(double dt) {
		int i,j,ij;
		double *ep=e;

		// Note: current structure of this loop is not suited to OpenMP.
		// It would need some changes to the variables for OpenMP.
		for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) {
			sum_e=(eta_(i+0.5,j)+eta_(i-0.5,j)+eta(i,j-0.5)+eta(i,j+0.5));
			ep[0]=eta_(i+0.5,j)
			ep[1]=eta_(i-0.5,j)
			ep[2]=eta_(i,j-0.5)
			ep[3]=eta_(i,j+0.5)
			ep[4]=sum_e;
			ep[5]=1/sum_e;
			ep+=6;
		}
	}
	/** Function to determine whether a grid point is on the edge or not.
	 */
	inline bool edge(int i,int ij) {return i==0||i==m-1||ij>=mn-m||ij<m;}
	/** Functions to specify the corner stencil entries. */
	inline double a_dl(int i,int ij) {return 0;}
	inline double a_dr(int i,int ij) {return 0;}
	inline double a_ul(int i,int ij) {return 0;}
	inline double a_ur(int i,int ij) {return 0;}
	/** Functions to specify the vertical stencil entries. */
	inline double a_dc(int i,int ij) {return edge(i,ij)?ep[2]:0;}
	inline double a_uc(int i,int ij) {return edge(i,ij)?ep[3]:0;}
	/** Functions to specify the horizontal stencil entries. */
	inline double a_cl(int i,int ij) {return edge(i,ij)?ep[1]:0;}
	inline double a_cr(int i,int ij) {return edge(i,ij)?ep[0]:0;}
	/** Function to specify the central stencil entry (on the diagonal of
	 * the linear system). */
	inline double a_cc(int i,int ij) {return ep[4];}
	/** Function to multiply by the reciprocal of the central stencil entry.
	 * This is specified as a separate function for computational
	 * efficiency. */
	inline double inv_cc(int i,int ij,double v) {return ep[5]*v;}
	/** Calculates the ith component of the multiplication (A-D)z, needed
	 * in the Gauss--Seidel smoothing iteration. */
	inline double mul_a(int i,int ij) {
		return edge(i,ij)?0:ep[0]*z[ij+1]+ep[1]*z[ij-1]+ep[3]*z[ij+m]+ep[2]*z[ij-m];

	//	This function should return exactly this:
	//	return a_dl(i,ij)*z[ij-m-1]+a_dc(i,ij)*z[ij-m]+a_dr(i,ij)*z[ij-m+1]
	//	      +a_cl(i,ij)*z[ij-1]+a_cr(i,ij)*z[ij+1]
	//	      +a_ul(i,ij)*z[ij+m-1]+a_uc(i,ij)*z[ij+m]+a_ur(i,ij)*z[ij+m+1];
	}
};

int main() {
	const int m=1025,n=1025,mn=m*n;
	const double ax=-8,bx=8,ay=ax,by=bx;
	int i,j,ij;
	double *b=new double[mn],*z=new double[mn],x,y;
	multisetup1 msu(m,n,ax,bx,ay,by,z);
	tgmg<multisetup1,double,double> mg(msu,b,z);

	// Set up the solution and source arrays
	for(ij=j=0,y=ay;j<n;j++,y+=msu.dy) {
		for(i=0,x=ax;i<m;i++,x+=msu.dx,ij++) {
			z[ij]=0;
			b[ij]=msu.edge(i,ij)?0:1/(x*x+y*y+4);
		}
	}

	// Set up linear system coeffecients
	msu.create_table();

	// Set up the multigrid hierarchy
	mg.setup();

	// Carry out a multigrid solve
	mg.solve_gauss_seidel();

	// Output the solutions in a format that can be read by Gnuplot using
	// the command "splot 'filename' matrix binary"
	mg.output_b("b.0");
	mg.output_z("z.0");
	mg.output_res("r.0");

	// Delete dynamically allocated memory
	delete [] z;
	delete [] b;
}
