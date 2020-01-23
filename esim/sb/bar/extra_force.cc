#include "common.hh"
#include "vec.hh"
#include "extra_force.hh"

#include <cmath>

extra_force::extra_force(bd_sim &bds_) : bds(bds_) {
	bds.eforce=this;
	char buf[256];
	sprintf(buf,"%s/ffile",bds.filename);
	fp=safe_fopen(buf,"wb");
}

extra_force::~extra_force() {
	fclose(fp);
}

vec extra_force1::force(int i,int j) {
	double &t=bds.time;
	const double fr=0.25,frsq=fr*fr,fri=4*atan(1.0)/fr;
	double dx=bds.ax+i*bds.dx-x,dy=bds.ay+j*bds.dy-y,
	       rsq=dx*dx+dy*dy,fm=-1.85*5e-6*(t<4e5?t:(t>8e5?0:8e5-t)),//>0.45?0.5:bds.time),
	       o=rsq>frsq?0:-0.5*fm*(1+cos(sqrt(rsq)*fri));
	return vec(12*o,o);
}

vec extra_force2::force(int i,int j) {
	double &t=bds.time;
	const double fr=0.15,frsq=fr*fr,fri=4*atan(1.0)/fr;
	double dx=bds.ax+i*bds.dx-x,dy=bds.ay+j*bds.dy-y,
	       rsq=dx*dx+dy*dy,fm=-5*5e-6*t,o=rsq>frsq?0:-0.5*fm*(1+cos(sqrt(rsq)*fri));
	return vec(o,0);
}

vec extra_force3::force(int i,int j) {
	double &t=bds.time;
	const double fr=0.5,frsq=fr*fr,fri=4*atan(1.0)/fr;
	double dx=bds.ax+i*bds.dx-x,dy=bds.ay+j*bds.dy-y,
	       rsq=dx*dx+dy*dy,fm=-1.85*5e-6*(t<4e5?t:(t>8e5?0:8e5-t)),//>0.45?0.5:bds.time),
	       o=rsq>frsq?0:-0.5*fm*(1+cos(sqrt(rsq)*fri))*0.25;
	return vec(12*o,o);
}
