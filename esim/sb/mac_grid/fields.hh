#ifndef SHEAR_FIELDS_HH
#define SHEAR_FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct field {
	/** The horizontal velocity. */
	double u;
	/** The vertical velocity. */
	double v;
	/** A component of stress. */
	double s11;
	/** A component of stress. */
	double s12;
	/** A component of stress. */
	double s21;
	/** A component of stress. */
	double s22;
	/** A component of stress. */
	double s33;
	/** The x coordinate of the reference map. */
	double X;
	/** The y coordinate of the reference map. */
	double Y;
	/** The change in the horizontal velocity. */
	double cu;
	/** The change in the vertical velocity. */
	double cv;
	/** The change in a component of stress. */
	double cs11;
	/** The change in a component of stress. */
	double cs12;
	/** The change in a component of stress. */
	double cs21;
	/** The change in a component of stress. */
	double cs22;
	/** The change in a component of stress. */
	double cs33;
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
		s11+=cs11;
		s12+=cs12;
		s21+=cs21;
		s22+=cs22;
		s33+=cs33;
		X+=cX;
		Y+=cY;
	}
	inline void clear_stress() {
		s11=s12=s21=s22=s33=0;
	}
	inline void copy(field &f,double delx,double dely) {
		u=f.u;v=f.v;
		s11=f.s11;s12=f.s12;
		s21=f.s21;s22=f.s22;
		s33=f.s33;
		X=f.X+delx;Y=f.Y+dely;
	}
};

#endif
