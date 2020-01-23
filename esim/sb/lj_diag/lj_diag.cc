#include <cstdio>
#include <cstdlib>

#include "lj_sim.hh"

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i,j;
	srand(1729);

	lj_sim lj(512,5,5,5,0.1);
	lj.gravity=0.08;
	for(i=0;i<110;i++) lj.put(10*rnd()-5,60*rnd()-30,0.425+0.15*rnd());

	char buf[512];
	lj.compute_acceleration();
	//lj.output_gnuplot("test/g.0");
	for(i=1;i<=105;i++) {
		for(j=0;j<100;j++) lj.integrate(0.01);
	//	sprintf(buf,"test/g.%d",i);
	//	lj.output_gnuplot(buf);
	}

	lj.damp=0.2;
	lj.gravity=0;
	lj.label_atoms();
	lj.vy_fix=+0.01;
	for(i=0;i<400;i++) lj.integrate(0.01);
	lj.damp=1e-6;
	lj.vy_fix=0;
	for(i=0;i<100;i++) lj.integrate(0.01);
	lj.output_gnuplot("jiggle/g.0");
	for(i=1;i<=800;i++) {
		for(j=0;j<100;j++) lj.integrate(0.002);
		sprintf(buf,"jiggle/g.%d",i);
		lj.output_gnuplot(buf);
	}

	lj.damp=0.1;
	lj.vx_fix=0.0025;
	lj.vy_fix=0;
	lj.output_gnuplot("shear/g.0");
	for(i=1;i<=1000;i++) {
		for(j=0;j<1000;j++) lj.integrate(0.002);
		sprintf(buf,"shear/g.%d",i);
		lj.output_gnuplot(buf);
	}

	lj.damp=0.1;
	lj.vx_fix=0.0025;
	lj.vy_fix=0;
	lj.output_gnuplot("shear2/g.0");
	for(i=1;i<=700;i++) {
		for(j=0;j<700;j++) lj.integrate(0.002);
		sprintf(buf,"shear2/g.%d",i);
		lj.output_gnuplot(buf);
	}
}
