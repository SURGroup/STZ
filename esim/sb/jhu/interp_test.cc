#include <cstdio>
#include <cstdlib>
#include <cmath>

class interp {
	public:
		const int m,n,mn;
		const double ax,bx;
		const double ay,by;
		const double dx,dy;
		const double xsp,ysp;
		double *u;
		double a[16];
		int ijc;
		interp(int m_,int n_,double ax_,double bx_,double ay_,double by_)
			: m(m_), n(n_), mn(m_*n_), ax(ax_), bx(bx_), ay(ay_), by(by_),
			dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)), xsp(1/dx), ysp(1/dy),
			u(new double[mn]), ijc(-1) {
				for(int j=0;j<n;j++) for(int i=0;i<m;i++) {
					u[i+j*m]=cos(i+2*j);
				}
				a[0]=a[1]=a[2]=a[3]=0;
				a[4]=a[5]=a[6]=a[7]=0;
				a[8]=a[9]=a[10]=a[11]=0;
				a[12]=a[13]=a[14]=a[15]=0;
			}
		~interp() {
			delete [] u;
		}
		double f(double x,double y) {
			x=(x-ax)*xsp;y=(y-ay)*ysp;
			int i=static_cast<int>(x),j=static_cast<int>(y),ij;
			if(i<0) i=0;else if(i>m-2) i=m-2;
			if(j<0) j=0;else if(j>n-2) j=n-2;
			ij=i+m*j;
			if(ijc!=ij) table_setup(i,j,ij);
			return f_int(x-i,y-j);
		}
		void table_setup(int i,int j,int ij) {
			ijc=ij;
			double *up=u+ij;
			double c00,c01,c02,c03;
			double c10,c11,c12,c13;
			double c20,c21,c22,c23;
			double c30,c31,c32,c33;
			compute_x(i,up,c01,c11,c21,c31);
			compute_x(i,up+m,c02,c12,c22,c32);
			if(j==0) {
				compute_x(i,up+2*m,c03,c13,c23,c33);
				fill_ad(a,c01,c02,c03);
				fill_ad(a+4,c11,c12,c13);
				fill_ad(a+8,c21,c22,c23);
				fill_ad(a+12,c31,c32,c33);
			} else if(j==n-2) {
				compute_x(i,up-m,c00,c10,c20,c30);
				fill_au(a,c00,c01,c02);
				fill_au(a+4,c10,c11,c12);
				fill_au(a+8,c20,c21,c22);
				fill_au(a+12,c30,c31,c32);
			} else {
				compute_x(i,up-m,c00,c10,c20,c30);
				compute_x(i,up+2*m,c03,c13,c23,c33);
				fill_a(a,c00,c01,c02,c03);
				fill_a(a+4,c10,c11,c12,c13);
				fill_a(a+8,c20,c21,c22,c23);
				fill_a(a+12,c30,c31,c32,c33);
			}
		}
		double f_int(double x,double y) {
			double z;
			z=x;x=y;y=z;
			return fx(a,x)+y*(fx(a+4,x)+y*(fx(a+8,x)+y*fx(a+12,x)));
		}
	private:
		inline double fx(double *ap,double y) {
			return *ap+y*(ap[1]+y*(ap[2]+y*ap[3]));
		}
		void fill_ad(double *ap,double c1,double c2,double c3) {
			*ap=c1;
			ap[1]=-1.5*c1+2*c2-0.5*c3;
			ap[2]=0.5*c1-c2+0.5*c3;
			ap[3]=0;
		}
		void fill_au(double *ap,double c0,double c1,double c2) {
			*ap=c1;
			ap[1]=-0.5*c0+0.5*c2;
			ap[2]=0.5*c0-c1+0.5*c2;
			ap[3]=0;
		}
		void fill_a(double *ap,double c0,double c1,double c2,double c3) {
			*ap=c1;
			ap[1]=-0.5*c0+0.5*c2;
			ap[2]=c0-2.5*c1+2*c2-0.5*c3;
			ap[3]=-0.5*c0+1.5*c1-1.5*c2+0.5*c3;
		}
		void compute_x(int i,double *up,double &c0,double &c1,double &c2,double &c3) {
			c0=*up;
			if(i==0) {
				c1=-1.5*(*up)+2*up[1]-0.5*up[2];
				c2=0.5*(*up)-up[1]+0.5*up[2];
				c3=0;
			} else if(i==m-2) {
				c1=-0.5*up[-1]+0.5*up[1];
				c2=0.5*up[-1]-*up+0.5*up[1];
				c3=0;
			} else {
				c1=-0.5*up[-1]+0.5*up[1];
				c2=up[-1]-2.5*(*up)+2*up[1]-0.5*up[2];
				c3=-0.5*up[-1]+1.5*(*up)-1.5*up[1]+0.5*up[2];
			}
		}
};

int main() {
	interp in(5,5,0,1,0,1);

	for(int j=0;j<29;j++) {
		double y=0.2+j*0.05;
		for(int i=0;i<29;i++) {
			double x=-0.2+i*0.05;
			printf("%g %g %g\n",x,y,in.f(x,y));
		}
	}
}
