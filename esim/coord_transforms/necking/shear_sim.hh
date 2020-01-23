#ifndef STZSIM_HH
#define STZSIM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

#include "fields.hh"
#include "stz_model.hh"
#include "qs_multi.hh"
#include "tgmg.hh"
#include "bi_interp.hh"

/** The number of standard deviations to consider when doing Gaussian
 * smoothing. */
const double gaussian_fac=5.;
using std::string;

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
        /** A time scale associated with the plastic deformation. */
        const double t_scale;
        /** A measure of what one unit of simulation time corresponds to in
         * seconds. This is used to properly scale the plastic deformation. */
        const double adapt_fac;
        /* The velocity with which we pull on the right boundary. */
        double v0;
        /** A pointer to the class to carry out the STZ plasticty. */
        stz_dynamics *stz;
        /** The filename of the output directory. */
        const char *filename;
        /** An array containing the simulation fields. */
        c_field* const fbase;
        /** The (0,0) position within the fbase array. */
        c_field* const fm;
        /** The current simulation time. */
        double time;
        /** The value of T_Z, used to convert the dimensionless effective
         * temperatures into Kelvins. */
        double TZ;
        /** The limiting value of effective temperature. */
        double chi_inf;
        /** The current output frame number. */
        int f_num;

        int curr_step=0;
        /* Pointer to an array holding the boundary terms and ghost regions. */
        bdry_field *bdrys_g;
        /* Pointer to the first physical boundary value in the boundary array. */
        bdry_field *bdrys;
        shear_sim(const int m_, const int n_, const double ax_, const double bx_,
              const double ay_, const double by_, const double mu_, const double K_,
              const double visc_, const double t_scale_, const double adapt_fac_, const double v0_,
              stz_dynamics *stz_, const unsigned int fflags_, const char *filename_); 
        ~shear_sim();
        void solve(double duration,int frames);
        void solve_quasistatic(double duration,int frames,int steps);
        void step_forward(double dt);
        void step_forward_quasistatic(double dt);
        void advection_step(double dt);
        void projection_step(double dt);
        void init_fields(int chi_case,double tem_base=600,double tem_delta=200);
        void write_files(int k);
        void initialize_random(double l,double tem_base=600,double tem_delta=200);
        double gaussian_normal_factor(double llinv,int cut);

        /* Transform functions. */
        inline mat3 calc_Ts(double L, double dC_dX, double Y, double dW_dX, double W);
        inline mat3 calc_Ts_inv_trans(double L, double dC_dX, double Y, double dW_dX, double W);
        inline mat3 calc_Ts_dot(double Ldot, double dC_dX_dot, double Y, double dW_dX_dot, double VY, double dW_dX, double Wdot);
        inline mat3 calc_l_trans(c_field &f, c_field &fr, c_field &fu, c_field &fur, 
                                 bdry_field &bf, bdry_field &bfr, double L, double Ldot, double X, double Y,
                                 mat3 &Ts_inv_trans, int i, int j);
        inline mat3 calc_sig(mat3 &Ts, mat3 &Sig) { return Ts*Sig*Ts.transpose(); };
        inline mat3 calc_sig0(mat3 &sig) { return sig - mat3(sig.trace()/3.); };
        inline double calc_sbar(mat3 &sig0);
        
        // Pulling on the boundary.
        //inline double calc_L()     { return (time < 1)? 1 + .5*v0*time*time : 1 + v0*time; }
        //inline double calc_Ldot()  { return (time < 1)? v0*time             : v0;      }
        //inline double calc_Lddot() { return (time < 1)? v0                  : 0;       }

        inline double calc_L()     { return 1 + v0*time; }
        inline double calc_Ldot()  { return v0;      }
        inline double calc_Lddot() { return 0;       }

    private:
        /** An object containing the configuration of the linear system
         * to be solved using the multigrid method. */
        qs_multi qsm;
        void set_boundaries();
        void output(const char *prefix,const int mode,const int sn);
        inline void net_force(double &fx,double &fy,double xf,double yf,c_field *fp);
        inline void rm_eno2(double &Xd,double &Yd,double hs,c_field *fp,int d);
        inline void bdry_adv_eno2(double &Wadv_d, double &Cadv_d, double hs, bdry_field *fp, int d);
        inline void top_bdry_eno2(double &p_d, double &p_dd, double hs, bdry_field *fp, int d);
        inline void bot_bdry_eno2(double &p_d, double &p_dd, double hs, bdry_field *fp, int d);
        inline void rmv_eno2(double &ud,double &vd,double &Xd,double &Yd,double hs,c_field *fp,int d);
        inline st_field st_eno2(double hs,c_field *fp,int d);
        inline double eno2(double p0,double p1,double p2,double p3);
        void update_bdry_profiles(double dt, double hx);
        double adaptive_plastic_term(double sbar, double &chiv, double &ddev, double dt);
        inline int step_mod(int a,int b) {return a>=0?a%b:b-1-(b-1-a)%b;}
        void box_muller(double &r0, double &r1);
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
