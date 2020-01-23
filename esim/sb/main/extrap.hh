#ifndef SHEAR_EXTRAP_HH
#define SHEAR_EXTRAP_HH

#include "fields.hh"
#include "level++.hh"

class ex_ref_map {
	public:
		const int c;
		ex_ref_map(c_field *fm_) : c(2), fm(fm_) {}
		inline void acc(int ij,double *p,double w,double x,double y) {
			c_field &f=fm[ij];
			*p+=w*f.X;p[1]+=w*x*f.X;p[2]+=w*y*f.X;
			p[3]+=w*f.Y;p[4]+=w*x*f.Y;p[5]+=w*y*f.Y;
		}
		inline void set(int ij,double *p,double al,double be,double ga) {
			c_field &f=fm[ij];
			f.X=*p*al+p[1]*be+p[2]*ga;
			f.Y=p[3]*al+p[4]*be+p[5]*ga;
		}
	private:
		c_field* fm;
};

class ex_staggered {
	public:
		const int c;
		ex_staggered(c_field *fm_) : c(4), fm(fm_) {}
		inline void acc(int ij,double *p,double w,double x,double y) {
			c_field &f=fm[ij];
			*p+=w*f.u;p[1]+=w*x*f.u;p[2]+=w*y*f.u;
			p[3]+=w*f.v;p[4]+=w*x*f.v;p[5]+=w*y*f.v;
			p[6]+=w*f.X;p[7]+=w*x*f.X;p[8]+=w*y*f.X;
			p[9]+=w*f.Y;p[10]+=w*x*f.Y;p[11]+=w*y*f.Y;
		}
		inline void set(int ij,double *p,double al,double be,double ga) {
			c_field &f=fm[ij];
			f.u=*p*al+p[1]*be+p[2]*ga;
			f.v=p[3]*al+p[4]*be+p[5]*ga;
			f.X=p[6]*al+p[7]*be+p[8]*ga;
			f.Y=p[9]*al+p[10]*be+p[11]*ga;
		}
	private:
		c_field* fm;
};

class ex_regular {
	public:
		const int c;
		ex_regular(c_field *fm_,int m_) : c(5), fm(fm_), m(m_) {}
		inline void acc(int ij,double *p,double w,double x,double y) {
			c_field &f=fm[r(ij)];
			*p+=w*f.p;p[1]+=w*x*f.p;p[2]+=w*y*f.p;
			p[3]+=w*f.q;p[4]+=w*x*f.q;p[5]+=w*y*f.q;
			p[6]+=w*f.s;p[7]+=w*x*f.s;p[8]+=w*y*f.s;
			p[9]+=w*f.tau;p[10]+=w*x*f.tau;p[11]+=w*y*f.tau;
			p[12]+=w*f.chi;p[13]+=w*x*f.chi;p[14]+=w*y*f.chi;
		}
		inline void set(int ij,double *p,double al,double be,double ga) {
			c_field &f=fm[r(ij)];
			f.p=*p*al+p[1]*be+p[2]*ga;
			f.q=p[3]*al+p[4]*be+p[5]*ga;
			f.s=p[6]*al+p[7]*be+p[8]*ga;
			f.tau=p[9]*al+p[10]*be+p[11]*ga;
			f.chi=p[12]*al+p[13]*be+p[14]*ga;
		}
	private:
		c_field* fm;
		const int m;
		inline int r(int ij) {return ij+ij/m;}
};

#endif
