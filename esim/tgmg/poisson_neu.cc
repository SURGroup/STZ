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
	const double fex,fey;
	/** A pointer to the solution vector. */
	double* const z;
	multisetup1(const int m_,const int n_,const double ax_,const double bx_,const double ay_,const double by_,double* const z_)
		: m(m_), n(n_), mn(m_*n_), x_prd(false), y_prd(false),
		acc(1e-20), ax(ax_), bx(bx_), ay(ay_), by(by_),
		dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), dxdy(dx/dy), dydx(dy/dx),
		fex(1/(dx*dx)), fey(1/(dy*dy)), z(z_) {}
	/** Functions to specify the corner stencil entries. */
	inline double a_dl(int i,int ij) {return 0;}
	inline double a_dr(int i,int ij) {return 0;}
	inline double a_ul(int i,int ij) {return 0;}
	inline double a_ur(int i,int ij) {return 0;}
	/** Functions to specify the vertical stencil entries. */
	inline double a_dc(int i,int ij) {return ij<m?0:fey;}
	inline double a_uc(int i,int ij) {return ij>=mn-m?0:fey;}
	/** Functions to specify the horizontal stencil entries. */
	inline double a_cl(int i,int ij) {return i==0?0:fex;}
	inline double a_cr(int i,int ij) {return i==m-1?0:fex;}
	/** Function to specify the central stencil entry (on the diagonal of
	 * the linear system). */
	inline double a_cc(int i,int ij) {return -a_dc(i,ij)-a_uc(i,ij)-a_cl(i,ij)-a_cr(i,ij);}
	/** Function to multiply by the reciprocal of the central stencil
	 * entry. This is specified as a separate function for computational
	 * efficiency. */
	inline double inv_cc(int i,int ij,double v) {return v/a_cc(i,ij);}
	/** Calculates the ith component of the multiplication (A-D)z, needed
	 * in the Gauss--Seidel smoothing iteration. */
	inline double mul_a(int i,int ij) {
		return fey*((ij<m?0:z[ij-m])+(ij>=mn-m?0:z[ij+m]))
		      +fex*((i==0?0:z[ij-1])+(i==m-1?0:z[ij+1]));
	}
};

int main() {
	const int m=129,n=129,mn=m*n;
	const double ax=-8,bx=8,ay=ax,by=bx;
	int i,j,ij;
	double *b=new double[mn],*z=new double[mn],x,y;
	multisetup1 msu(m,n,ax,bx,ay,by,z);
	tgmg<multisetup1,double,double> mg(msu,b,z);

	// Set up the solution and source arrays
	for(ij=j=0,y=ay;j<n;j++,y+=msu.dy) {
		for(i=0,x=ax;i<m;i++,x+=msu.dx,ij++) {
			z[ij]=0;
			b[ij]=0;
		}
	}

	b[20+m*20]=1;
	b[5+m*5]=-1;

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
