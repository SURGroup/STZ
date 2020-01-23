#ifndef C3D_MAT3_HH
#define C3D_MAT3_HH

#include <cstdio>
#include <cmath>
#include <limits>
#include "vec3.hh"

class sym_mat3;

struct mat3 {
	double a11,a12,a13;
	double a21,a22,a23;
	double a31,a32,a33;
	mat3() {};
	mat3(double a) : a11(a), a12(0), a13(0), a21(0), a22(a), a23(0), a31(0), a32(0), a33(a) {};
	mat3(double a,double b,double c) : a11(a), a12(0), a13(0), a21(0), a22(b), a23(0), a31(0), a32(0), a33(c) {};
	mat3(double b11,double b12,double b13,double b21,double b22,double b23,double b31,double b32,double b33)
		: a11(b11), a12(b12), a13(b13), a21(b21), a22(b22), a23(b23), a31(b31), a32(b32), a33(b33) {};
	inline void set(double b11,double b12,double b13,double b21,double b22,double b23,double b31,double b32,double b33) {
		a11=b11;a12=b12;a13=b13;
		a21=b21;a22=b22;a23=b23;
		a31=b31;a32=b32;a33=b33;
	}

    // Returns an element according to indices (like a 2d index lookup).
    // Useful for looping, as necessary in 3d transformation code.
    double operator() (int ii, int jj) {
        switch(jj + 3*ii) {
            case 0:
                return a11;
                break;
            case 1:
                return a12;
                break;
            case 2:
                return a13;
                break;
            case 3:
                return a21;
                break;
            case 4:
                return a22;
                break;
            case 5:
                return a23;
                break;
            case 6:
                return a31;
                break;
            case 7:
                return a32;
                break;
            case 8:
                return a33;
                break;
            default:
                fprintf(stderr, "Invalid element access in mat3.ind()!!!\n");
                return std::numeric_limits<double>::quiet_NaN();
                break;
        }
    }

    // Sets an element according to indices.
    void ind_set(int ii, int jj, double val) {
        switch(jj + 3*ii) {
            case 0:
                a11 = val;
                break;
            case 1:
                a12 = val;
                break;
            case 2:
                a13 = val;
                break;
            case 3:
                a21 = val;
                break;
            case 4:
                a22 = val;
                break;
            case 5:
                a23 = val;
                break;
            case 6:
                a31 = val;
                break;
            case 7:
                a32 = val;
                break;
            case 8:
                a33 = val;
                break;
            default:
                fprintf(stderr, "Invalid element access in mat3.ind_set()!!!\n");
                *this = std::numeric_limits<double>::quiet_NaN();
                break;
        }
    }


    inline mat3 operator= (double val){
        a11 = a12 = a13 = a21 = a22 = a23 = a31 = a32 = a33 = val;
        return *this;
    }

	inline mat3 operator+ (mat3 p) {
		return mat3(a11+p.a11,a12+p.a12,a13+p.a13,
			    a21+p.a21,a22+p.a22,a23+p.a23,
			    a31+p.a31,a32+p.a32,a33+p.a33);
	}

	inline mat3 operator- (mat3 p) {
		return mat3(a11-p.a11,a12-p.a12,a13-p.a13,
			    a21-p.a21,a22-p.a22,a23-p.a23,
			    a31-p.a31,a32-p.a32,a33-p.a33);
	}

    inline void operator+= (mat3 p) {
        a11 += p.a11; a12 += p.a12; a13 += p.a13;
        a21 += p.a21; a22 += p.a22; a23 += p.a23;
        a31 += p.a31; a32 += p.a32; a33 += p.a33;
    }

    inline void operator-= (mat3 p) {
        a11 -= p.a11; a12 -= p.a12; a13 -= p.a13;
        a21 -= p.a21; a22 -= p.a22; a23 -= p.a23;
        a31 -= p.a31; a32 -= p.a32; a33 -= p.a33;
    }

    /* in-place Matrix-Scalar Multiplication. */
	inline void operator*= (double p) {
		a11*=p;a12*=p;a13*=p;
		a21*=p;a22*=p;a23*=p;
		a31*=p;a32*=p;a33*=p;
	}

    /* Matrix-Scalar Multiplication. */
	inline mat3 operator* (double p) {
		return mat3(a11*p,a12*p,a13*p,
			    a21*p,a22*p,a23*p,
			    a31*p,a32*p,a33*p);
	}

    /* Matrix-Vector Multiplication. */
    inline vec3 operator* (vec3 e) {
        return vec3(a11*e.x + a12*e.y + a13*e.z, 
                    a21*e.x + a22*e.y + a23*e.z, 
                    a31*e.x + a32*e.y + a33*e.z);
    }

    /* Matrix-Scalar Division. */
	inline mat3 operator/ (double p) {
		double pinv(1/p);
		return mat3(a11*pinv,a12*pinv,a13*pinv,
			    a21*pinv,a22*pinv,a23*pinv,
			    a31*pinv,a32*pinv,a33*pinv);
	}

    /* Matrix-Matrix Multiplication. */
	inline mat3 operator* (mat3 p) {
		return mat3(a11*p.a11+a12*p.a21+a13*p.a31,a11*p.a12+a12*p.a22+a13*p.a32,a11*p.a13+a12*p.a23+a13*p.a33,
			    a21*p.a11+a22*p.a21+a23*p.a31,a21*p.a12+a22*p.a22+a23*p.a32,a21*p.a13+a22*p.a23+a23*p.a33,
			    a31*p.a11+a32*p.a21+a33*p.a31,a31*p.a12+a32*p.a22+a33*p.a32,a31*p.a13+a32*p.a23+a33*p.a33);
	}

    /* Matrix-Matrix Multiplication. */
	inline mat3 operator* (sym_mat3 p);
	inline double det() {
		return a11*(a22*a33-a23*a32)+a21*(a13*a32-a12*a33)+a31*(a12*a23-a13*a22);
	}

	inline mat3 transpose() {
		return mat3(a11,a21,a31,a12,a22,a32,a13,a23,a33);
	}

	inline mat3 inverse(double &dt) {
		dt=det();
		double idet(1/dt);
		return mat3((a22*a33-a23*a32)*idet,(a13*a32-a12*a33)*idet,(a12*a23-a13*a22)*idet,
			    (a23*a31-a21*a33)*idet,(a11*a33-a13*a31)*idet,(a13*a21-a11*a23)*idet,
			    (a21*a32-a22*a31)*idet,(a12*a31-a11*a32)*idet,(a11*a22-a12*a21)*idet);
	}

	inline double trace() {
		return a11+a22+a33;
	}

	inline void print() {
		printf("%g %g %g\n%g %g %g\n%g %g %g\n",a11,a12,a13,a21,a22,a23,a31,a32,a33);
	}

	inline double modsq() {
		return 0.5*(a11*a11+a12*a12+a13*a13+a21*a21+a22*a22+a23*a23+a31*a31+a32*a32+a33*a33);
	}

	inline double mod() {
		return sqrt(modsq());
	}

	sym_mat3 trans_mult();
};

/* Scalar-First Scalar-Matrix Multiplication */
inline mat3 operator*(const double e, mat3 f) {
    return f*e;
}

/* Scalar-MatInv Multiplication */
inline mat3 operator/(double e, mat3 f) {
    double dt;
	return e*f.inverse(dt);
}

class sym_mat3 {
	public:
		double a11, a12, a13;
		double a22, a23;
		double a33;

		sym_mat3() {};

		sym_mat3(double a_): a11(a_), a12(0), a13(0), a22(a_), a23(0), a33(a_) {};

        sym_mat3(double a, double b, double c) : a11(a), a12(0), a13(0), a22(b), a23(0), a33(c) {};

		sym_mat3(double b11, double b12, double b13, double b22, double b23, double b33)
			: a11(b11), a12(b12), a13(b13), 
              a22(b22), a23(b23), a33(b33) {};

        sym_mat3(double b11, double b12, double b13, 
                 double b21, double b22, double b23, 
                 double b31, double b32, double b33)
            : a11(b11), a12(b12), a13(b13), 
              a22(b22), a23(b23), a33(b33) {};

		inline void put(double *p) {
			*(p++) = a11; *(p++) = a12; *(p++) = a13;
            *(p++) = a22; *(p++) = a23; *p = a33;
		}

        inline sym_mat3 operator= (double val){
            a11 = a12 = a13 = a22 = a23 = a33 = val;
            return *this;
        }

		inline sym_mat3 operator+ (sym_mat3 p) {
			return sym_mat3(a11 + p.a11, a12 + p.a12, a13 + p.a13,
                            a22 + p.a22, a23 + p.a23, a33 + p.a33);
		}

        inline mat3 operator+ (mat3 p){
            return mat3(a11 + p.a11, a12 + p.a12, a13 + p.a13,
                        a12 + p.a21, a22 + p.a22, a23 + p.a23,
                        a13 + p.a31, a23 + p.a32, a33 + p.a33);
        }

        inline void operator+= (sym_mat3 p){
            a11 += p.a11; a12 += p.a12; a13 += p.a13;
            a22 += p.a22; a23 += p.a23; a33 += p.a33;
        }

		inline sym_mat3 operator- (sym_mat3 p) {
			return sym_mat3(a11 - p.a11, a12 - p.a12, a13 - p.a13,
                            a22 - p.a22, a23 - p.a23, a33 - p.a33);
		}

        inline mat3 operator- (mat3 p){
            return mat3(a11 - p.a11, a12 - p.a12, a13 - p.a13,
                        a12 - p.a21, a22 - p.a22, a23 - p.a23,
                        a13 - p.a31, a23 - p.a32, a33 - p.a33);
        }

        inline void operator-= (sym_mat3 p){
            a11 -= p.a11; a12 -= p.a12; a13 -= p.a13;
            a22 -= p.a22; a23 -= p.a23; a33 -= p.a33;
        }

		inline sym_mat3 operator* (double p) {
			return sym_mat3(a11*p, a12*p, a13*p,
                            a22*p, a23*p, a33*p);
		}

        inline void operator*= (double p){
            a11 *= p; a12 *= p; a13 *= p;
            a22 *= p; a23 *= p; a33 *= p;
        }

        /* Symmetric Matrix-Vector Multiplication */
        inline vec3 operator* (vec3 e) {
        return vec3(a11*e.x + a12*e.y + a13*e.z, 
                    a12*e.x + a22*e.y + a23*e.z, 
                    a13*e.x + a23*e.y + a33*e.z);
        }

		inline sym_mat3 operator/ (double p) {
			double pinv(1./p);
			return sym_mat3(a11*pinv, a12*pinv, a13*pinv,
                            a22*pinv, a23*pinv, a33*pinv);
		}

        inline sym_mat3 inverse(double &dt) {
		dt = det();
		double idet(1./dt);
		return sym_mat3(
                (a22*a33-a23*a23)*idet,(a13*a23-a12*a33)*idet,(a12*a23-a13*a22)*idet,
			    (a23*a13-a12*a33)*idet,(a11*a33-a13*a13)*idet,(a13*a12-a11*a23)*idet,
			    (a12*a23-a22*a13)*idet,(a12*a13-a11*a23)*idet,(a11*a22-a12*a12)*idet);
        }

		inline double trace() {
			return a11 + a22 + a33;
		}

		inline double det() {
			return a11*(a22*a33 - a23*a23) + a12*(2*a13*a23 - a12*a33) - a13*a13*a22;
		}

        inline sym_mat3 transpose() {
            return *this;
        }

		sym_mat3 trans_mult();

		inline void pack(double *&p) {
			*(p++) = a11; *(p++) = a12; *(p++) = a13;
			*(p++) = a22; *(p++) = a23; *(p++) = a33;
		}

		inline void print() {
			printf("%g %g %g\n%g %g %g\n%g %g %g\n", a11, a12, a13, a12, a22, a23, a13, a23, a33);
		}

		inline void clear() {
			a11 = a12 = a13 = a22 = a23 = a33 = 0;
		}

        inline double modsq() {
            return 0.5*(a11*a11 + a12*a12 + a13*a13 + a12*a12 + a22*a22 + a23*a23 + a13*a13 + a23*a23 + a33*a33);
        }

        inline double mod() {
            return sqrt(modsq());
        }

        mat3 operator*(mat3 p){
		return mat3(
                a11*p.a11 + a12*p.a21 + a13*p.a31, a11*p.a12 + a12*p.a22 + a13*p.a32, a11*p.a13 + a12*p.a23 + a13*p.a33,
			    a12*p.a11 + a22*p.a21 + a23*p.a31, a12*p.a12 + a22*p.a22 + a23*p.a32, a12*p.a13 + a22*p.a23 + a23*p.a33,
			    a13*p.a11 + a23*p.a21 + a33*p.a31, a13*p.a12 + a23*p.a22 + a33*p.a32, a13*p.a13 + a23*p.a23 + a33*p.a33);
        }

        mat3 operator*(sym_mat3 p){
		return mat3(
                a11*p.a11 + a12*p.a12 + a13*p.a13, a11*p.a12 + a12*p.a22 + a13*p.a23, a11*p.a13 + a12*p.a23 + a13*p.a33,
			    a12*p.a11 + a22*p.a12 + a23*p.a13, a12*p.a12 + a22*p.a22 + a23*p.a23, a12*p.a13 + a22*p.a23 + a23*p.a33,
			    a13*p.a11 + a23*p.a12 + a33*p.a13, a13*p.a12 + a23*p.a22 + a33*p.a23, a13*p.a13 + a23*p.a23 + a33*p.a33);
        }

		void project_to_other(double v1,double v2,double v3,sym_mat3 &s);
		void project(double v1,double v2,double v3);
		void eigenvalues(double &eig1,double &eig2,double &eig3);
		void eigenvectors(double &eig1,double &eig2,double &eig3,mat3 &Lam);

	private:
		void compute_eigenvector(double eig,double &a,double &b,double &c,double t);
		inline bool nonzero(double a,double threshold) {
			return std::abs(a)>threshold;
		}
};

/*
 *inline mat3 operator*(const double p, mat3 m) {
 *    return mat3(p*m.a11,p*m.a12,p*m.a13,
 *            p*m.a21,p*m.a22,p*m.a23,
 *            p*m.a31,p*m.a32,p*m.a33);
 *}
 */


/* Double-First Double-sym_mat3 Multiplication */
inline sym_mat3 operator*(const double p, sym_mat3 m) {
	return sym_mat3(p*m.a11, p*m.a12, p*m.a13,
                    p*m.a22, p*m.a23, p*m.a33);
}

/* Scalar-MatInv Multiplication */
inline sym_mat3 operator/(double e, sym_mat3 f) {
    double dt;
	return e*f.inverse(dt);
}

inline mat3 mat3::operator* (sym_mat3 p) {
	return mat3(a11*p.a11+a12*p.a12+a13*p.a13,a11*p.a12+a12*p.a22+a13*p.a23,a11*p.a13+a12*p.a23+a13*p.a33,
		    a21*p.a11+a22*p.a12+a23*p.a13,a21*p.a12+a22*p.a22+a23*p.a23,a21*p.a13+a22*p.a23+a23*p.a33,
		    a31*p.a11+a32*p.a12+a33*p.a13,a31*p.a12+a32*p.a22+a33*p.a23,a31*p.a13+a32*p.a23+a33*p.a33);
}

#endif
