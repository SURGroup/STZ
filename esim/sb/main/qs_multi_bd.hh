#ifndef QS_MULTI_BD_HH
#define QS_MULTI_BD_HH

#include "vec.hh"
#include "mat.hh"
#include "fields.hh"

class bd_sim;
class extra_force;

class qs_multi_bd {
	public:
		const int m,ml,n,mn;
		const bool x_prd;
		const bool y_prd;
		static const char gs_mode=0;
		c_field* const fm;
		vec* const z;
		vec* const b;
		extra_force *&eforce;
		double *phi;
		int *c;
		int *cc;
		const double acc;
		const double K,qK,Kmu,mu;
		const double dx,dy,xsp,ysp,xsp2,ysp2,xyfac;
		const double tmufac1,tmufac2;
		const double bcs;
		const double ocs;
		const double ncs;
		double idt,Kmx,Kmy,mx,my,cxx,cyy;
		const mat rz,rc,mrc;
		mat rh,rv,ru,rb,ro,rn,iu,ib,io,in;
		int *d;
		qs_multi_bd(bd_sim &bds);
		~qs_multi_bd();
		void init(double visc,double dt);
		inline mat a_dl(int i,int ij) {return d[ij]==-1?rc:(d[ij]>=0?a[d[ij]]:rz);}
		inline mat a_dc(int i,int ij) {return d[ij]==-1?rv:(d[ij]>=0?a[d[ij]+1]:rz);}
		inline mat a_dr(int i,int ij) {return d[ij]==-1?mrc:(d[ij]>=0?a[d[ij]+2]:rz);}
		inline mat a_cl(int i,int ij) {return d[ij]==-1?rh:(d[ij]>=0?a[d[ij]+3]:rz);}
		inline mat a_cc(int i,int ij) {return d[ij]==-1?ru:(d[ij]==-3?ro:(d[ij]==-2?rb:rn));}
		inline mat a_cr(int i,int ij) {return d[ij]==-1?rh:(d[ij]>=0?a[d[ij]+4]:rz);}
		inline mat a_ul(int i,int ij) {return d[ij]==-1?mrc:(d[ij]>=0?a[d[ij]+5]:rz);}
		inline mat a_uc(int i,int ij) {return d[ij]==-1?rv:(d[ij]>=0?a[d[ij]+6]:rz);}
		inline mat a_ur(int i,int ij) {return d[ij]==-1?rc:(d[ij]>=0?a[d[ij]+7]:rz);}
		inline vec inv_cc(int i,int ij,vec v) {return (d[ij]==-1?iu:(d[ij]==-3?io:(d[ij]==-2?ib:in)))*v;}
		inline vec mul_a(int i,int ij) {
			if(d[ij]==-1) {
				vec* zp=z+ij,*zd=zp-1,*zu=zp+1;
				return vec(Kmx*(zd->x+zu->x)+my*(zp[-m].x+zp[m].x)
				   +qK*(zd[-m].y+zu[m].y-zu[-m].y-zd[m].y),
				   mx*(zd->y+zu->y)+Kmy*(zp[-m].y+zp[m].y)
				   +qK*(zd[-m].x+zu[m].x-zu[-m].x-zd[m].x));
			}
			if(d[ij]<0) return vec(0);
			mat *ap=a+d[ij];
			vec *zp=z+ij;
			return *ap*zp[-m-1]+ap[1]*zp[-m]+ap[2]*zp[-m+1]+ap[3]*zp[-1]
				+ap[4]*zp[1]+ap[5]*zp[m-1]+ap[6]*zp[m]+ap[7]*zp[m+1];
		}
	private:
		void square1(double z,mat &h,mat &v,mat &a0,mat &a1,mat &a2,mat &a3,bool pos);
		void square2(double z,mat &h,mat &v,mat &a0,mat &a1,mat &a2,mat &a3,bool pos);
		inline double f(double zz) {return zz<0.5?0:(zz>1.5?2:2*(zz-0.5));}
		vec contrib(int i,int j,double nx,double ny);
		int na,a_mem;
		mat *a;
		void add_normals_memory();
};

#endif
