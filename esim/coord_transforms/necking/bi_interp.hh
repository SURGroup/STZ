#ifndef BI_INTERP_HH
#define BI_INTERP_HH

#include "fields.hh"

class shear_sim;

class bicubic_interp {
    public:
        shear_sim &ss;
        const bool cor;
        const int mode;
        const int m,n;
        const double ax,ay;
        const double xsp,ysp;
        double a[16];
        int ijc;
        bicubic_interp(shear_sim &ss_,int mode_);
        ~bicubic_interp() {}
        double f(double x,double y);
        double f_grad_f(double x,double y,double &fx,double &fy);
    private:
        void table_setup(int i,int j,int ij);
        void compute_x(int i,c_field *up,double &c0,double &c1,double &c2,double &c3);
        inline double yl(double *ap,double y) {
            return *ap+y*(ap[1]+y*(ap[2]+y*ap[3]));
        }
        inline double dyl(double *ap,double y) {
            return ap[1]+y*(2*ap[2]+3*y*ap[3]);
        }
        void grid_index(double &x,double &y);
        void fill_a(double *ap,double c0,double c1,double c2,double c3);
};

#endif
