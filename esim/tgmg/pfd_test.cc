#include <cmath>

#include "tgmg.hh"
#include "poisson_fd.hh"

int main() {
    const int m=1025,n=1025;
    const double ax=-8,bx=8,ay=ax,by=bx;
    poisson_fd pfd(m,n,ax,bx,ay,by);

    // Initialize the multigrid solver, set the verbosity to maximum, and set
    // up the multigrid hierarchy
    tgmg<poisson_fd,double,double> mg(pfd,pfd.b,pfd.z);
    mg.verbose=3;
    mg.setup();

    // Set up the solution and source arrays
    mg.clear_z();
    pfd.gaussian_source_term(0,0,1,1);

    // Solve using multigrid V-cycles
    mg.solve_v_cycle();

    // Output the solutions in a format that can be read by Gnuplot using
    // the command "splot 'filename' matrix binary"
    const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);
    mg.output_b("b.0",ax,dx,ay,dy);
    mg.output_z("z.0",ax,dx,ay,dy);
}
