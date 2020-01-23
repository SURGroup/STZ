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
	/** Computes the magnitude squared of the deviatoric stress
	 * tensor.
	 * \param[in] ij the gridpoint to consider.
	 * \return The magnitude squared. */
	inline double devsq() {
		return (3*q*q+s*s+tau*tau);
	}
	/** Computes the magnitude of the deviatoric stress tensor.
	 * \param[in] ij the gridpoint to consider.
	 * \return The magnitude. */
	inline double dev() {
		return sqrt(devsq());
	}
	inline double sxx() {
		return -p+s-q;
	}
	inline double syy() {
		return -p-s-q;
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
