#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#define SLIDING_BOUNDARY
#define IMPLICIT
#define LINEAR_EXTRAPOLATION
//#define FULL_ENO
#define VISCOSITY 0.002

const double horiz_speed=0.005;
const double vert_speed=0;

#include "vec.hh"
#include "../levelset/level++.cc"
#include "sbsim.cc"

const double bw=0.1;
const double bd=0.15;
const double init_w_pos=0.4;

inline double max(double a,double b) {return a>b?a:b;}

int main() {
//	sbsim sim(300,150,-1,1,-149/299.0,149.0/299.0,0.15,1,1e8,30,60,1,"output9");
	sbsim sim(300,400,-299.0/399.0,299.0/399.0,-299/399.0,299.0/399.0,0.15,1,1e9,30,60,2,"fat4");
//	sbsim sim(300,150,-1,1,-199/399.0,199.0/399.0,0.15,1,1e10,30,60,2,"output35");
//	sbsim sim(500,250,-1,1,-249/499.0,249.0/499.0,0.15,1,1e8,30,60,1,"output15");
//	sim.wallx=0.6;sim.wallu=horiz_speed;sim.wallv=0;
	sim.wallx=init_w_pos;sim.wallu=horiz_speed;sim.wallv=vert_speed;
	sim.solve(0,50,500);
}

void sbsim::init_fields() {
	int i,j,ij;double x,y,yy;
	for(ij=j=0;j<n;j++) {
		y=ay+dy*j;yy=abs(y)-0.500001;
		for(i=0;i<m;i++,ij++) {
			x=ax+dx*i;
			phi[ij]=abs(x)<bw?max(yy,y-(0.500001-bd*(0.5+0.5*cos(x*pi/bw)))):yy;
			v[ij]=0;//abs(y)<0.400001?(x>(init_w_pos+1e-5)?vert_speed:(x<-(init_w_pos+1e-5)?-vert_speed:0)):0;
			u[ij]=0;
			p[ij]=0;
			s[ij]=tau[ij]=0;
			chi[ij]=0.074;
		}
	}
}
