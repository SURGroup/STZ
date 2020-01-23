#ifndef SHEAR_FIELDS_HH
#define SHEAR_FIELDS_HH

#include <cmath>
#include "mat3.hh"
using std::isnan;

/** Data structure for storing the fields at grid points. */
struct c_field {
    /** The horizontal velocity. */
    double u;
    /** The vertical velocity. */
    double v;
    /* 11 component of stress. */
    double s11;
    /* 12 component of stress. */
    double s12;
    /* 22 component of stress. */
    double s22;
    /* 33 component of stress. */
    double s33;
    /** The effective temperature. */
    double chi;
    /** The x coordinate of the reference map. */
    double X;
    /** The y coordinate of the reference map. */
    double Y;
    /** The adaptively-computed D_pl term. */
    double ddev;
    /* sbar. */
    double dev;
	/** The change in the horizontal velocity. (Note that the same variable
	 * is also used to saving the total deformation rate while extracting
	 * different output at the end of each frame step.) */
	double cu;
	/** The change in the vertical velocity. (Note that the same variable
	 * is also used to save the total deformation rate while extracting
	 * different output at the end of each frame step.) */
	double cv;
    /* The change in the 11 component of stress. */
	double cs11;
    /* The change in the 12 component of stress. */
	double cs12;
    /* The change in the 22 component of stress. */
	double cs22;
    /* The change in the 33 component of stress. */
	double cs33;
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
        u     += cu;
        v     += cv;
        cu = cv = 0;
        update_others();
    }
    inline void update_staggered() {
        s11 += cs11;
        s12 += cs12;
        s22 += cs22;
        s33 += cs33;
        chi = cchi;
        cs11 = cs12 = cs22 = cs33 = cchi = 0;
    }
    inline void update_ref_map() {
        X += cX;
        Y += cY;
    }
    inline void update_regular() {
        u += cu;
        v += cv;
        update_ref_map();
    }
    inline void update_others() {
        update_ref_map();
        update_staggered();
    }
    inline void extrap(c_field &f0, c_field &f1) {
        u   = 2*f0.u   - f1.u;
        v   = 2*f0.v   - f1.v;
        s11 = 2*f0.s11 - f1.s11;
        s12 = 2*f0.s12 - f1.s12;
        s22 = 2*f0.s22 - f1.s22;
        s33 = 2*f0.s33 - f1.s33;
        chi = 2*f0.chi - f1.chi;
    }
    inline mat3 Sigma() {
        return mat3(s11, s12,   0,
                    s12, s22,   0,
                      0,   0, s33);
    }
    inline double fval(int i) {
        switch(i) {
            case 0: return u;
            case 1: return v;
            case 2: return s11;
            case 3: return s12;
            case 4: return s22;
            case 5: return s33;
            case 6: return chi;
            case 7: return X;
            case 8: return Y;
            case 9: return cu;
            case 10: return cv;
            case 11: return cs11;
            case 12: return cs12;
            case 13: return cs22;
            default: return cs33;
        }
    }
    /** Diagnostic function for finding weird values. */
    inline bool weird() {
        bool rslt = false;
        if(isnan(u)||fabs(u)>1e10)     { printf("u is weird!\n");   rslt = true; }
        if(isnan(v)||fabs(v)>1e10)     { printf("v is weird!\n");   rslt = true; }
        if(isnan(s11)||fabs(s11)>1e10) { printf("s11 is weird!\n"); rslt =  true; }
        if(isnan(s12)||fabs(s12)>1e10) { printf("s12 is weird!\n"); rslt =  true; }
        if(isnan(s22)||fabs(s22)>1e10) { printf("s22 is weird!\n"); rslt =  true; }
        if(isnan(s33)||fabs(s33)>1e10) { printf("s33 is weird!\n"); rslt =  true; }
        if(isnan(chi)||fabs(chi)>1e10) { printf("chi is weird!\n"); rslt =  true; }
        return rslt;
    }
};

/** Data structure for storing the fields at grid points. */
struct st_field {
    /** s11 component of stress. */
    double s11;
    /** s12 component of stress. */
    double s12;
    /** s22 component of stress. */
    double s22;
    /** s33 component of stress. */
    double s33;
    /** The effective temperature. */
    double chi;
    st_field() {}
    st_field(double s11_, double s12_, double s22_, double s33_, double chi_)
        : s11(s11_), s12(s12_), s22(s22_), s33(s33_), chi(chi_) {}
};

/** Data structure for storing the boundary profiles and their advective derivatives. */
struct bdry_field {
    /* Top profile */
    double tp;
    /* Top profile advective time derivative. */
    double tp_adv;
    /* Bottom profile */
    double bp;
    /* Bottom profile advective time derivative. */
    double bp_adv;
    /* Change in top profile. */
    double ctp;
    /* Change in top profile advective time derivative. */
    double ctp_adv;
    /* Change in bottom profile. */
    double cbp;
    /* Change in bottom profile advective time derivative. */
    double cbp_adv;
    /* Merge the updates for the boundary profiles. */
    inline void update() {
        tp     += ctp;
        tp_adv += ctp_adv;
        bp     += cbp;
        bp_adv += cbp_adv;
        
        ctp = ctp_adv = cbp = cbp_adv = 0;
    }
    // Compute W from T and B.
    inline double comp_W() {
        //double tmp = xx*xx - 1;
        //return 1 - tt*tmp*tmp;
        return .5*(tp - bp);
    }
    // Compute the advective derivative of W from the advective derivatives of T and B.
    inline double comp_Wadv() {
        //double tmp = xx*xx - 1;
        //return -tmp*tmp - 4*xx*tt*tmp*u;
        return .5*(tp_adv - bp_adv);
    }
    // Compute C from T and B.
    inline double comp_C() {
        //return 0;
        return .5*(tp + bp);
    }
    // Compute the advective derivative of C from the advective derivatives of T and B.
    inline double comp_Cadv() {
        //return 0;
        return .5*(tp_adv + bp_adv);
    }
    // Linear extrapolation.
    inline void extrap(bdry_field &f0, bdry_field &f1) {
        tp     = 2*f0.tp     - f1.tp;
        tp_adv = 2*f0.tp_adv - f1.tp_adv;
        bp     = 2*f0.bp     - f1.bp;
        bp_adv = 2*f0.bp_adv - f1.bp_adv;
    }
    bdry_field() : tp(1), tp_adv(0), bp(-1), bp_adv(0), ctp(0), ctp_adv(0), cbp(0), cbp_adv(0) {}
};

#endif
