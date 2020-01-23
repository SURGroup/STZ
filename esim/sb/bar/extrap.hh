#ifndef SHEAR_EXTRAP_HH
#define SHEAR_EXTRAP_HH

#include "fields.hh"
#include "level++.hh"

class ex_ref_map {
	public:
		const int c;
		ex_ref_map(c_field *fm_) : c(2), fm(fm_) {}
		inline void dset(double *&fp,double fac,int ijl,int ijr) {
			c_field &fl=fm[ijl],&fr=fm[ijr];
			*(fp++)=fac*(fr.X-fl.X);*(fp++)=fac*(fr.Y-fl.Y);
		}
		inline void set1st(double *&fp,int ij) {
			c_field &f=fm[ij];
			*(fp++)=f.X;*(fp++)=f.Y;
		}
		inline void set2nd(double *&fp,int ij1,int ij2) {
			c_field &f1=fm[ij1],&f2=fm[ij2];
			*(fp++)=2*f1.X-0.5*f2.X;
			*(fp++)=2*f1.Y-0.5*f2.Y;
		}
		inline void set(int ij,double *a0,double b0,double *a1,double b1,
				double *a2,double b2,double *a3,double b3) {
			c_field &f=fm[ij];
			f.X=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.Y=*a0*b0+*a1*b1+*a2*b2+*a3*b3;
		}
	private:
		c_field* fm;
};

class ex_staggered {
	public:
		const int c;
		ex_staggered(c_field *fm_) : c(4), fm(fm_) {}
		inline void dset(double *&fp,double fac,int ijl,int ijr) {
			c_field &fl=fm[ijl],&fr=fm[ijr];
			*(fp++)=fac*(fr.u-fl.u);*(fp++)=fac*(fr.v-fl.v);
			*(fp++)=fac*(fr.X-fl.X);*(fp++)=fac*(fr.Y-fl.Y);
		}
		inline void set1st(double *&fp,int ij) {
			c_field &f=fm[ij];
			*(fp++)=f.u;*(fp++)=f.v;
			*(fp++)=f.X;*(fp++)=f.Y;
		}
		inline void set2nd(double *&fp,int ij1,int ij2) {
			c_field &f1=fm[ij1],&f2=fm[ij2];
			*(fp++)=2*f1.u-0.5*f2.u;
			*(fp++)=2*f1.v-0.5*f2.v;
			*(fp++)=2*f1.X-0.5*f2.X;
			*(fp++)=2*f1.Y-0.5*f2.Y;
		}
		inline void set(int ij,double *a0,double b0,double *a1,double b1,
				double *a2,double b2,double *a3,double b3) {
			c_field &f=fm[ij];
			f.u=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.v=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.X=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.Y=*a0*b0+*a1*b1+*a2*b2+*a3*b3;
		}
	private:
		c_field* fm;

};

class ex_regular {
	public:
		const int c;
		ex_regular(c_field *fm_,int m_) : c(5), fm(fm_), m(m_) {}
		inline void dset(double *&fp,double fac,int ijl,int ijr) {
			c_field &fl=fm[r(ijl)],&fr=fm[r(ijr)];
			*(fp++)=fac*(fr.p-fl.p);*(fp++)=fac*(fr.q-fl.q);
			*(fp++)=fac*(fr.s-fl.s);*(fp++)=fac*(fr.tau-fl.tau);
			*(fp++)=fac*(fr.chi-fl.chi);
		}
		inline void set1st(double *&fp,int ij) {
			c_field &f=fm[r(ij)];
			*(fp++)=f.p;*(fp++)=f.q;*(fp++)=f.s;
			*(fp++)=f.tau;*(fp++)=f.chi;
		}
		inline void set2nd(double *&fp,int ij1,int ij2) {
			c_field &f1=fm[r(ij1)],&f2=fm[r(ij2)];
			*(fp++)=2*f1.p-0.5*f2.p;
			*(fp++)=2*f1.q-0.5*f2.q;
			*(fp++)=2*f1.s-0.5*f2.s;
			*(fp++)=2*f1.tau-0.5*f2.tau;
			*(fp++)=2*f1.chi-0.5*f2.chi;
		}
		inline void set(int ij,double *a0,double b0,double *a1,double b1,
				double *a2,double b2,double *a3,double b3) {
			c_field &f=fm[r(ij)];
			f.p=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.q=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.s=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.tau=*(a0++)*b0+*(a1++)*b1+*(a2++)*b2+*(a3++)*b3;
			f.chi=*a0*b0+*a1*b1+*a2*b2+*a3*b3;
		}
	private:
		c_field* fm;
		const int m;
		inline int r(int ij) {return ij+ij/m;}
};

#endif
