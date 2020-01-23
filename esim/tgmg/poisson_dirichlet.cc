// Copyright 2015 Chen-Hung Wu and Chris Rycroft,
// This is an example file for testing the multigrid code with Dirichlet B.C.

#include "tgmg.cc"

// Multisetup structure for Poisson problem
struct multisetup {
  const int m;
  const int n;
  const int mn;
  const bool x_prd;
  const bool y_prd;
  static const bool gs_mode=0;
  const double acc;
  const double ax, bx;
  const double ay, by;
  const double dx, dy;
  const double dxdy, dydx;
  const double fm, fm_inv, fex, fey, fc;
  double* const z;
  multisetup(const int m_, const int n_,
              const double ax_, const double bx_,
              const double ay_, const double by_, double* const z_)
    : m(m_), n(n_), mn(m_*n_), x_prd(false), y_prd(false),
    acc(1e-20), ax(ax_), bx(bx_), ay(ay_), by(by_),
    dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), dxdy(dx/dy), dydx(dy/dx),
    fm(-4/(dx*dx)), fm_inv(-0.25*dx*dx), fex(1/(dx*dx)),
    fey(1/(dy*dy)), fc(0), z(z_) {}

  inline bool edge(int i, int ij) {
    return i == m-1 || ij >= mn-m;
  }
  inline double a_dl(int i, int ij) {return 0;}
  inline double a_dr(int i, int ij) {return 0;}
  inline double a_ul(int i, int ij) {return 0;}
  inline double a_ur(int i, int ij) {return 0;}
  inline double a_dc(int i, int ij) {return edge(i,ij)?0:(ij<m?0:fey);}
  inline double a_uc(int i, int ij) {return edge(i,ij)?0:(ij>=mn-m?0:fey);}
  inline double a_cl(int i, int ij) {return edge(i,ij)?0:(i==0?0:fex);}
  inline double a_cr(int i, int ij) {return edge(i,ij)?0:(i==m-1?0:fex);}
  inline double a_cc(int i, int ij) {
	  return edge(i,ij)?-2*(fex+fey):-a_dc(i,ij)-a_uc(i,ij)-a_cl(i,ij)-a_cr(i,ij);
  }
  inline double inv_cc(int i, int ij, double v) {return v/a_cc(i,ij);}
  inline double mul_a(int i, int ij) {
    if(edge(i,ij)) return 0;
    return fey*((ij<m?0:z[ij-m])+(ij>=mn-m?0:z[ij+m]))
          +fex*((i==0?0:z[ij-1])+(i==m-1?0:z[ij+1]));
  }
};

int main() {
  const int m = 33, n = 33, mn = m*n;
  const double ax = -8, bx = 8, ay = ax, by = bx;
  int i, j, ij;
  double* b = new double[mn];
  double* z = new double[mn];
  double x, y;
  multisetup msu(m, n, ax, bx, ay, by, z);
  tgmg<multisetup, double, double> mg(msu, b, z);

  // Set up the solution and source arrays
  for (ij = j = 0, y = ay; j < n; j++, y += msu.dy) {
    for (i = 0, x = ax; i < m; i++, x += msu.dx, ij++) {
      z[ij] = 0;
      //b[ij] = msu.edge(i, ij) ? 0 : 1/(x*x+y*y+4);

      b[ij] = std::abs(x) > 4.0 || std::abs(y) > 4.0 ? 0 : -1;
    }
  }

  // Set up the multigrid hierarchy
  mg.setup();

  // Carry out ten multigrid solves for timing purposes
  mg.solve_v_cycle();

  // Optional output routines
  mg.output_b("b.0");
  mg.output_z("z.0");
  mg.output_res("r.0");

  // Delete dynamically allocated memory
  delete [] z;
  delete [] b;
}
