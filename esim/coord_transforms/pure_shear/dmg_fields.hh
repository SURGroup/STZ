#ifndef DMG_FIELDS_HH
#define DMG_FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct c_field_dmg {
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
	/** The grain size field. */
	double aa;
	/** The x coordinate of the reference map. */
	double X;
	/** The y coordinate of the reference map. */
	double Y;
	/** The change in the horizontal velocity. */
	double cu;
	/** The change in the vertical velocity. */
	double cv;
	/** The change in pressure. */
	double cp;
	/** The change in a component of shear stress. */
	double cs;
	/** The change in a component of shear stress. */
	double cq;
	/** The change in a component of shear stress. */
	double ctau;
	/** The change in effective temperature. */
	double cchi;
	/** The change in the grain size field. */
	double caa;
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
		aa=caa;
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
};

/** Data structure for storing the fields at grid points. */
struct st_field_dmg {
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
	double aa;
	st_field_dmg() {}
	st_field_dmg(double p_,double q_,double s_,double tau_,double chi_,double aa_)
		: p(p_), q(q_), s(s_), tau(tau_), chi(chi_), aa(aa_) {}
};

#endif
