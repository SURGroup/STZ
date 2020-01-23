#ifndef EXTRA_FORCE_CC
#define EXTRA_FORCE_CC

#include "bd_sim.hh"

class extra_force {
	public:
		double x;
		double y;
		extra_force(bd_sim &bds_);
		virtual ~extra_force();
		void write(int k) {
			fprintf(fp,"%d %g %g %g\n",k,bds.time,x,y);
			fflush(fp);
		}
		inline void compute(double dt) {
			bds.bicubic_velocity(x,y,cx,cy);
			cx*=dt;cy*=dt;
		}
		inline void update() {
			x+=cx;y+=cy;
		}
		virtual vec force(int i,int j) = 0;
	protected:
		bd_sim &bds;
		FILE *fp;
		double cx;
		double cy;
};

class extra_force1 : public extra_force {
	public:
		extra_force1(bd_sim &bds_) : extra_force(bds_) {x=-1;y=0;}
		virtual vec force(int i,int j);
};

class extra_force2 : public extra_force {
	public:
		extra_force2(bd_sim &bds_) : extra_force(bds_) {x=y=0;}
		virtual vec force(int i,int j);
};

class extra_force3 : public extra_force {
	public:
		extra_force3(bd_sim &bds_) : extra_force(bds_) {x=-1;y=0;}
		virtual vec force(int i,int j);
};

#endif
