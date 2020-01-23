#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "../shared/vec.hh"
#include "level++.hh"

const double pi=3.1415926535897932384626433832795;

inline double fphi(double x,double y) {return sqrt(x*x+y*y)-2;}
inline double fu(double x,double y) {return 2*x-y;}//-sin(3*x)+cos(3*y);}
inline double fd(double x,double y) {return x-2*y;}//-sin(3*x)+cos(3*y);}

int main(int argc,char **argv) {
	if(argc!=2) {
		fputs("One argument required\n",stderr);
		return 1;
	}
	int n=atof(argv[1]),ne(n+1),nne(ne*ne);
	double a=-(1-1.0/n)*pi,b=-a;
	double d=(b-a)/(n-1);

	levelset ls(n,n,a,a,d,d,-5,5);

	double x,y,*phi=ls.phi;
	int i,j,ij,k,l,*s=ls.s;

	// Set up level set function
	for(ij=j=0;j<n;j++) {
		y=a+d*j;
		for(i=0;i<n;i++,ij++) {
			x=a+d*i;
			phi[ij]=fphi(x,y);
		}
	}
	output("phi",0,ls);
	ls.build_band();

	// Set up staggered field within the body
	double *u=new double[nne],*up=u;
	ex_single e(u);

	int *ssp=ls.ss;
	for(j=0;j<ne;j++) {
		y=-pi+d*j;
		for(i=0;i<ne;i++,up++,ssp++) {
			ij=i+n*j;
			k=l=0;
			if(i>0) {
				if(j>0) {k++;if(s[ij-n-1]>3) l++;}
				if(j<n) {k++;if(s[ij-1]>3) l++;}
			}
			if(i<n) {
				if(j>0) {k++;if(s[ij-n]>3) l++;}
				if(j<n) {k++;if(s[ij]>3) l++;}
			}
			if(l==0) {
				*up=fu(-pi+d*i,y);
				*ssp=0;
			} else if(l==k) {
				*up=fd(-pi+d*i,y);
				*ssp=7;
			} else {
				*ssp=3;
				*up=0;
			}

		}
	}

	output("u",0,u,ne,ne,-pi,-pi,d,d);
	output("sg",0,ls.sg,ne,ne,-pi,-pi,d,d);
	output("ss",0,ls.ss,ne,ne,-pi,-pi,d,d);

//	ls.extrapolate_staggered_fields_reverse(e);
	ls.extrapolate_staggered_fields(e);

	output("uf",0,u,ne,ne,-pi,-pi,d,d);
	output("sgf",0,ls.sg,ne,ne,-pi,-pi,d,d);

	delete [] u;
}
