#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <limits>
#include <cmath>

const double pi=3.1415926535897932384626433832795;

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);

/** \brief Calculates the argument of two-dimensional position vector.
 *
 * Calculates the argument of the two-dimensional position vector.
 * \param[in] (x,y) the coordinates of the vector.
 * \return The argument. */
inline double argument(double x,double y) {
	const double pi=3.1415926535897932384626433832795;
	return x+y>0?(x>y?atan(y/x):0.5*pi-atan(x/y)):(x>y?-atan(x/y)-0.5*pi:atan(y/x)+(y>0?pi:-pi));
}

#endif
