#ifndef SHEAR_FIELDS_HH
#define SHEAR_FIELDS_HH

#include <cmath>
using std::isnan;

/** Data structure for storing the fields at grid points. */
struct c_field {
    /** The horizontal velocity. */
    double u;
    /** The vertical velocity. */
    double v;
    /** The pressure. */
    double p;
    /** A component of shear stress. */
    double q;
    /** A component of shear stress. */
    double s;
    /** A component of shear stress. */
    double tau;
    /** The effective temperature. */
    double chi;
    /** The x coordinate of the reference map. */
    double X;
    /** The y coordinate of the reference map. */
    double Y;
    /** The adaptively-computed D_pl term. */
    double ddev;
    /** The change in the horizontal velocity. (Note that the same variable
     * is also used to saving the total deformation rate while extracting
     * different output at the end of each frame step.) */
    double cu;
    /** The change in the vertical velocity. (Note that the same variable
     * is also used to saving the total deformation rate while extracting
     * different output at the end of each frame step.) */
    double cv;
    /** The change in pressure. (Note that the same variable is also used
     * to saving the total deformation rate while extracting different
     * output at the end of each frame step.) */
    double cp;
    /** The change in a component of shear stress. */
    double cs;
    /** The change in a component of shear stress. */
    double cq;
    /** The change in a component of shear stress. */
    double ctau;
    /** The change in effective temperature. */
    double cchi;
    /** The change in the x coordinate of the reference map. */
    double cX;
    /** The change in the y coordinate of the reference map. */
    double cY;
    /** Applies the stored updates to the fields. For velocities, pressure,
     * and the reference map, the updates are added. For the shear stresses
     * and effective temperature (which are semi-implicitly updated), the
     * values are replaced with the stored updates. */
    inline void update() {
        u+=cu;
        v+=cv;
        update_others();
    }
    inline void update_staggered() {
        p+=cp;
        q+=cq;
        s+=cs;
        tau+=ctau;
        chi=cchi;
    }
    inline void update_ref_map() {
        X+=cX;
        Y+=cY;
    }
    inline void update_regular() {
        u+=cu;
        v+=cv;
        update_ref_map();
    }
    inline void update_others() {
        update_ref_map();
        update_staggered();
    }
    inline void prd_bc(c_field &f,double dX,double dY) {
        u=f.u;v=f.v;p=f.p;
        q=f.q;s=f.s;tau=f.tau;
        chi=f.chi;X=f.X+dX;Y=f.Y+dY;
    }
    inline void set_corner(double u_,double v_,double X_,double Y_) {
        u=u_;v=v_;X=X_;Y=Y_;
    }
    inline void extrap(c_field &f0,c_field &f1) {
        p=2*f0.p-f1.p;q=2*f0.q-f1.q;
        s=2*f0.s-f1.s;tau=2*f0.tau-f1.tau;
        chi=2*f0.chi-f1.chi;
    }
    inline void strain_extrap(c_field &f0,c_field &f1) {
        cu=2*f0.cu-f1.cu;cv=2*f0.cv-f1.cv;
        cp=2*f0.cp-f1.cp;cq=2*f0.cq-f1.cq;
        cs=2*f0.cs-f1.cs;ctau=2*f0.ctau-f1.ctau;
    }
    inline void strain_prd_bc(c_field &f) {
        cu=f.cu;cv=f.cv;cp=f.cp;
        cq=f.cq;cs=f.cs;ctau=f.ctau;
    }
    /** Computes the magnitude squared of the deviatoric stress tensor, taking
     * into account the transformation of the grid.
     * \param[in] lt the simple shear transformation amount.
     * \return The magnitude squared. */
    inline double devsq(double lt) {

/*        Transformed Calculation
 *
 *
 *        double lt  = lt;
 *        double lt2 = lt*lt;
 *        double lt3 = lt*lt2;
 *        double lt4 = lt2*lt2;
 *
 *        double reg      = 3*q*q + s*s + tau*tau;
 *        double lt_term  = 2*lt*(p + 2*(q + s))*tau;
 *        double lt2_term = lt2*third*(3*p*p + 6*q*q + 6*q*s + 3*p*(3*q + s) + 4*tau*tau);
 *        double lt3_term = lt3*2*third*(p + q + s)*tau;
 *        double lt4_term = lt4*third*(p + q + s)*(p + q + s);
 *
 *        return reg + lt_term + lt2_term + lt3_term + lt4_term;
 *
 */

        // First compute the untransformed coordinates.
        double unt_tau = tau - lt*(p + q + s);
        double unt_q   = q + lt*(1/3.)*(-tau + lt*.5*(p + q + s));
        double unt_s   = s + lt*(tau - lt*.5*(p + q + s));

        // Then calculate sbar using the untransformed coordinates.
        return 3*unt_q*unt_q + unt_s*unt_s + unt_tau*unt_tau;

        /* Old calculation */
        //return (3*q*q+s*s+tau*tau);
    }
    /** Computes the magnitude of the deviatoric stress tensor, taking into
     * account the transformation of the grid.
     * \param[in] lt the simple shear transformation amount.
     * \return The magnitude. */
    inline double dev(double lt) {
        return sqrt(devsq(lt));
    }
    inline double sxx() {
        return -p+s-q;
    }
    inline double syy() {
        return -p-s-q;
    }
    inline double fval(int i) {
        switch(i) {
            case 0: return u;
            case 1: return v;
            case 2: return p;
            case 3: return q;
            case 4: return s;
            case 5: return tau;
            case 6: return chi;
            case 7: return X;
            case 8: return Y;
            case 9: return cu;
            case 10: return cv;
            case 11: return cp;
            case 12: return cq;
            case 13: return cs;
            default: return ctau;
        }
    }
    /** Diagnostic function for finding weird values. */
    inline bool weird() {
        if(isnan(u)||fabs(u)>1e10) return true;
        if(isnan(v)||fabs(v)>1e10) return true;
        if(isnan(p)||fabs(p)>1e10) return true;
        if(isnan(q)||fabs(q)>1e10) return true;
        if(isnan(s)||fabs(s)>1e10) return true;
        if(isnan(tau)||fabs(tau)>1e10) return true;
        if(isnan(chi)||fabs(chi)>1e10) return true;
        return false;
    }

};

/** Data structure for storing the fields at grid points. */
struct st_field {
    /** The pressure. */
    double p;
    /** A component of shear stress. */
    double q;
    /** A component of shear stress. */
    double s;
    /** A component of shear stress. */
    double tau;
    /** The effective temperature. */
    double chi;
    st_field() {}
    st_field(double p_,double q_,double s_,double tau_,double chi_)
        : p(p_), q(q_), s(s_), tau(tau_), chi(chi_) {}
};

#endif
