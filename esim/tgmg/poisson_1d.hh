#ifndef POISSON_1D_HH
#define POISSON_1D_HH

#include "conj_grad.hh"

class poisson_1d : public conj_grad {
    public:
        /** The grid spacing. */
        const double dx;
        /** The inverse grid spacing squared. */
        const double xsp2;
        poisson_1d(int dof_);
        void init();
        void print_solution();
        virtual void mul_A(double *in,double *out);
};

#endif
