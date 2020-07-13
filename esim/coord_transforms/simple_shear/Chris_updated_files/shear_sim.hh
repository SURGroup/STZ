#ifndef STZSIM_HH
#define STZSIM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "vec.hh"
#include "mat.hh"
#include "fields.hh"
#include "stz_model.hh"
#include "qs_multi.hh"
#include "tgmg.hh"
#include "bi_interp.hh"

/** The number of standard deviations to consider when doing Gaussian
 * smoothing. */
const double gaussian_fac=5.;

/** \brief A class to carry out a 2D elasticity simulation. */
class shear_sim {
    public:
        /** The number of interior grid points in the horizontal direction. */
        const int m;
        /** The number of grid points in the horizontal direction, including
         * ghost regions. */
        const int gm;
        /** The number of interior grid points in the vertical direction. */
        const int n;
        /** The number of grid points in the vertical direction in the linear
         * system. */
        const int lsn;
        /** The number of grid points in the vertical direction, including
         * ghost regions. */
        const int gn;
        /** The total number of interior grid points. */
        const int mn;
        /** The total number of grid points, including grid points */
        const int gmn;
        /** Output fields and miscelleneous flags. */
        const unsigned int fflags;
        /** The lower bound in the x direction. */
        const double ax;
        /** The lower bound in the y direction. */
        const double ay;
        /** The upper bound in the x direction. */
        const double bx;
        /** The upper bound in the y direction. */
        const double by;
        /** The width of the simulation. */
        const double lx;
        /** The height of the simulation. */
        const double ly;
        /** The grid spacing in the x direction. */
        const double dx;
        /** The grid spacing in the y direction. */
        const double dy;
        /** The inverse grid spacing in the x direction. */
        const double xsp;
        /** The inverse grid spacing in the y direction. */
        const double ysp;
        /** The elastic shear modulus. */
        const double mu;
        /** The reciprocal of the elastic shear modulus. */
        const double mu_inv;
        /** The elastic bulk modulus. */
        const double K;
        /** The viscosity. */
        const double viscosity;
        /** The amount of Dpl-mediated effective temperature diffusion. */
        const double chi_len;
        /** A time scale associated with the plastic deformation. */
        const double t_scale;
        /** A measure of what one unit of simulation time corresponds to in
         * seconds. This is used to properly scale the plastic deformation.
         */
        const double adapt_fac;
        /** The velocity of the top boundary. */
        const double u_bdry;
        /** The factor appearing in the exponents of A and B**/
        const double lamb;
        /** A pointer to the class to carry out the STZ plasticty. */
        stz_dynamics *stz;
        /** Whether or not the simulation is periodic in the y direction. */
        bool y_prd;
        /** The filename of the output directory. */
        const char *filename;
        /** An array containing the simulation fields. */
        c_field* const fbase;
        /** The (0,0) position within the fbase array. */
        c_field* const fm;
        /** An array for the Lagrangian tracer positions. */
        double* tr;
        /** The current simulation time. */
        double time;
        /** The value of T_Z, used to convert the dimensionless effective
         * temperatures into Kelvins. */
        double TZ;
        /** The limiting value of effective temperature. */
        double chi_inf;
        /** The current output frame number. */
        int f_num;
        /** The current number of tracers. */
        int ntrace;
        /** The size of the tracer grid. */
        int tr_m,tr_n;
        /** The time at which the tracers were last updated or created. */
        double tr_time;
        shear_sim(const int m_, const int n_, const double ax_, const double bx_,
              const double ay_, const double by_, const double mu_, const double K_,
              const double visc_, const double chi_len_, const double t_scale_,
              const double adapt_fac, const double u_bdry_, const double lamb_,
              stz_dynamics *stz_, bool y_prd_, const unsigned int fflags_,
              const char *filename_); ~shear_sim();
        void solve(double duration,int frames);
        void solve_quasistatic(double duration,int frames,int steps);
        void step_forward(double dt);
        void step_forward_quasistatic(double dt);
        void advection_step(double dt);
        void chi_diffusion(double dt);
        void projection_step(double dt);
        void init_fields(int chi_case,double tem_base=600,double tem_delta=200);
        void write_files(int k);
        void compute_strain();
        void tractions(double &bx,double &by,double &tx,double &ty);
        void l2_comparison(shear_sim &ss, double *l2);
        void l2_comparison_transform(shear_sim &ss, double *l2);
        void bilin_interp(double ii, double jj, double (&results)[7]);
        void initialize_chi_bicubic(int mm,int nn,double *arr,double lo=0,double sca=1);
        void initialize_random(double l,double tem_base=600,double tem_delta=200);
        void initialize_tracers(int mm,int nn);
        inline void initialize_tracers() {initialize_tracers(m,n);}
        void output_tracers(int sn);
        void output_tracers_matrix(const char *prefix,const int mode,const int sn);
        void update_tracers();
        double gaussian_normal_factor(double llinv,int cut);
        /** Returns the velocity of the boundary, based on the current
         * simulation time.
         * \return The velocity. */
        inline double bdry_vel() {
            //return time<1?time*u_bdry:u_bdry;
            return u_bdry;
        }
        /** Returns the displacement of the boundary, based on the
         * current simulation time.
         * \return The displacement. */
        inline double bdry_pos() {return time<1?u_bdry*time*time*0.5:u_bdry*(time-0.5);}
        inline double calc_A() {return exp(lamb*time);}
        inline double calc_dA_dt() {return lamb*exp(lamb*time);}
        inline double calc_dAsq_dtsq() {return lamb*lamb*exp(lamb*time);}
        inline double calc_B() {return exp(-lamb*time);}
        inline double calc_dB_dt() {return -lamb*exp(-lamb*time);}
        inline double calc_dBsq_dtsq() {return lamb*lamb*exp(-lamb*time);}
        inline double calc_dlogA_dt() {return lamb;}
        inline double calc_dlogB_dt() {return -lamb;}
        inline double calc_xi_x(int xx){ return lamb*(ax + xx*dx); }
        inline double calc_xi_y(int yy){ return -lamb*(ay + yy*dy); }

        /* Code to compute the untransformed stressed in terms of the transformed
         * stresses. See Eq. 84-87 in the document. */
        inline double calc_p(c_field &f) {
            return f.p + (1./3)*lamb*time*(-2*f.tau + lamb*time*(f.p + f.q + f.s));
        }
        /*
         *inline double calc_tau(c_field &f){
         *    return f.tau - lamb*time*(f.p + f.q + f.s);
         *}
         *inline double calc_q(c_field &f){
         *    double third = 1./3;
         *    return f.q - lamb*time*(f.tau*third + lamb*time*.5*(f.p + f.q + f.s));
         *}
         *inline double calc_s(c_field &f){
         *    return f.s + lamb*time*(f.tau - lamb*time*.5*(f.p + f.q + f.s));
         *}
         */
    private:
        /** An object containing the configuration of the linear system
         * to be solved using the multigrid method. */
        qs_multi qsm;
        /** Bicubic interpolation classes. */
        bicubic_interp *bi[10];
        void set_strain_boundaries();
        void set_boundaries();
        void output(const char *prefix,const int mode,const int sn);
        inline void net_force(double &fx,double &fy,double xf,double yf,c_field *fp);
        inline void rm_eno2(double &Xd,double &Yd,double hs,c_field *fp,int d);
        inline void rmv_eno2(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field *fp,int d);
        inline st_field st_eno2(double hs,c_field *fp,int d);
        inline double eno2(double p0,double p1,double p2,double p3);
        double adaptive_plastic_term(double sbar, double &chiv, double &ddev, double dt);
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        void box_muller(double &r0, double &r1);
        bool correct_tracer(double xref,double yref,double &xx,double &yy);
        void check_grid_wtf();
        /** Temporary storage for used during the output routine. */
        float *buf;
#ifdef _OPENMP
        inline double wtime() {return omp_get_wtime();}
#else
        inline double wtime() {return 0;}
#endif
};

#endif
