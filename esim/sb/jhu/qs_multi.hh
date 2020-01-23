#ifndef QS_MULTI_HH
#define QS_MULTI_HH

#include "vec.hh"
#include "mat.hh"
#include "tgmg.hh"

class shear_sim;

struct qs_multi {
	const int m,n,mn;
	const bool x_prd;
	const bool y_prd;
	static const char gs_mode=0;
	vec* const z;
	const double K,qK,Kmu,mu;
	const double dx,dy,xsp,ysp,xsp2,ysp2,xyfac;
	const double tmufac1,tmufac2;
	const double bcs;
	double viscdt_prev;
	double Kmx,Kmy,mx,my,cxx,cyy;
	const mat regz;
	const double acc;
	mat regh,regv,regc,regci,regcb,regcbi;
	qs_multi(shear_sim &ss);
	~qs_multi() {}
	void init(double viscdt);
	inline bool bdry(int i,int ij) {return ij>=mn-m||ij<m;}
	inline mat a_dl(int i,int ij) {return bdry(i,ij)?regz:mat(0,qK,qK,0);}
	inline mat a_dr(int i,int ij) {return bdry(i,ij)?regz:mat(0,-qK,-qK,0);}
	inline mat a_ul(int i,int ij) {return bdry(i,ij)?regz:mat(0,-qK,-qK,0);}
	inline mat a_ur(int i,int ij) {return bdry(i,ij)?regz:mat(0,qK,qK,0);}
	inline mat a_cl(int i,int ij) {return bdry(i,ij)?regz:regh;}
	inline mat a_cr(int i,int ij) {return bdry(i,ij)?regz:regh;}
	inline mat a_dc(int i,int ij) {return bdry(i,ij)?regz:regv;}
	inline mat a_uc(int i,int ij) {return bdry(i,ij)?regz:regv;}
	inline mat a_cc(int i,int ij) {return bdry(i,ij)?regcb:regc;}
	inline vec inv_cc(int i,int ij,vec v) {return (bdry(i,ij)?regcbi:regci)*v;}
	inline vec mul_a(int i,int ij) {
		if(bdry(i,ij)) return vec(0,0);
		vec* zp=z+ij,*zd=zp+(i==0?m-1:-1),*zu=zp+(i==m-1?1-m:1);
		return vec(Kmx*(zd->x+zu->x)+my*(zp[-m].x+zp[m].x)
			   +qK*(zd[-m].y+zu[m].y-zu[-m].y-zd[m].y),
			   mx*(zd->y+zu->y)+Kmy*(zp[-m].y+zp[m].y)
			   +qK*(zd[-m].x+zu[m].x-zu[-m].x-zd[m].x));
		//return a_dc(i,ij)*z[ij-m]+a_cl(i,ij)*z[ij-1]+a_cr(i,ij)*z[ij+1]+a_uc(i,ij)*z[ij+m]
		//	+a_dl(i,ij)*z[ij-m-1]+a_dr(i,ij)*z[ij-m+1]+a_ul(i,ij)*z[ij+m-1]+a_ur(i,ij)*z[ij+m+1];
	}
	inline void solve() {
		if(!mg.solve_v_cycle(tp)) exit(1);
	}
	/** The multigrid solver. */
	tgmg<qs_multi,vec,mat> mg;
	tgmg_predict tp;
};

#endif
