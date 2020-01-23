#include <cstdlib>

#include "lj_sim.hh"

void fatal_error(const char *p,int status) {
	fprintf(stderr,"%s\n",p);
	exit(status);
}

FILE* safe_fopen(const char *fn,const char *mode) {
	FILE *fp=fopen(fn,mode);
	if(fp==NULL) {
		fprintf(stderr,"Can't open file '%s'\n",fn);
		exit(1);
	}
	return fp;
}

lj_sim::lj_sim(int mem_,double sx_,double sy_,double cutoff_,double damp_) :
	mem(mem_), sx(sx_), sy(sy_), cutoffsq(cutoff_*cutoff_), damp(damp_),
	vx_fix(0), vy_fix(0),
	x(new double[2*mem]), v(new double[2*mem]), a(new double[2*mem]),
	rad(new double[mem]), ty(new int[mem]), n(0) {}

lj_sim::~lj_sim() {
	delete [] ty;
	delete [] rad;
	delete [] a;
	delete [] v;
	delete [] x;
}

void lj_sim::integrate(double dt) {
	int i;
	for(i=0;i<2*n;i++) v[i]+=0.5*dt*a[i];
	for(i=0;i<n;i++) {
		x[2*i]+=dt*v[2*i];
		if(x[2*i]>sx) x[2*i]-=2*sx;
		else if(x[2*i]<-sx) x[2*i]+=2*sx;
		x[2*i+1]+=dt*v[2*i+1];
	}
	compute_acceleration();
	for(i=0;i<2*n;i++) v[i]+=0.5*dt*a[i];
}

void lj_sim::compute_acceleration() {
	int i,j;
	double fx,fy;

	// Clear accelerations
	for(i=0;i<2*n;i++) a[i]=0;

	// Pairwise forces
	for(i=1;i<n;i++) for(j=0;j<i;j++) if(pair_force(i,j,fx,fy)) {
		a[2*i]-=fx;a[2*i+1]-=fy;
		a[2*j]+=fx;a[2*j+1]+=fy;
	}

	// Extra forces
	for(i=0;i<n;i++) {
		a[2*i]-=damp*v[2*i];
		a[2*i+1]-=damp*v[2*i+1]+gravity*x[2*i+1];
	}

	fix_atoms();
}

void lj_sim::label_atoms() {
	for(int i=0;i<n;i++) {
		if(x[2*i+1]>3.6) ty[i]=1;
		else if(x[2*i+1]<-3.6) ty[i]=2;
	}
}

void lj_sim::fix_atoms() {
	for(int i=0;i<n;i++) if(ty[i]>0) {
		if(ty[i]==1) {v[2*i]=vx_fix;v[2*i+1]=vy_fix;}
		else {v[2*i]=-vx_fix;v[2*i+1]=-vy_fix;}
		a[2*i]=a[2*i+1]=0;
	}
}

void lj_sim::put(double px,double py,double pr) {
	x[2*n]=px;x[2*n+1]=py;
	v[2*n]=0;v[2*n+1]=0;
	a[2*n]=0;a[2*n+1]=0;
	rad[n]=pr;
	ty[n++]=0;
}

bool lj_sim::pair_force(int i,int j,double &fx,double &fy) {
	double dx=x[2*i]-x[2*j],dy=x[2*i+1]-x[2*j+1],r2;
	if(dx>sx) dx-=2*sx;
	else if(dx<-sx) dx+=2*sx;
	r2=dx*dx+dy*dy;
	if(r2<cutoffsq) {
		if(r2<0.4) r2=0.4;
	//	printf("%g %g %g ",rad[i],rad[j],r2);
		r2=(rad[i]+rad[j])*(rad[i]+rad[j])/r2;
		double r6=r2*r2*r2;
		fx=0.05*dx*r2*r6*(-12*r6+6);
		fy=0.05*dy*r2*r6*(-12*r6+6);
	//	printf("%g %g\n",fx,fy);
		return true;
	}
	return false;
}

void lj_sim::save_state(FILE *fp) {
	fwrite(&n,sizeof(double),1,fp);
	fwrite(x,sizeof(double),2*n,fp);
	fwrite(v,sizeof(double),2*n,fp);
	fwrite(a,sizeof(double),2*n,fp);
	fwrite(rad,sizeof(double),n,fp);
}

void lj_sim::load_state(FILE *fp) {
	double nt;
	if(fread(&nt,sizeof(double),1,fp)!=1) fatal_error("Can't read header",2);
	int nn=int(nt);
	if(nn<0) fatal_error("Negative number of particles",2);
	if(nn>mem) fatal_error("Not enough memory allocated",2);
	n=nn;
	if(fread(x,sizeof(double),2*n,fp)!=(unsigned int) 2*n)
		fatal_error("Can't read position data",2);
	if(fread(v,sizeof(double),2*n,fp)!=(unsigned int) 2*n)
		fatal_error("Can't read velocity data",2);
	if(fread(a,sizeof(double),2*n,fp)!=(unsigned int) 2*n)
		fatal_error("Can't read acceleration data",2);
	if(fread(rad,sizeof(double),n,fp)!=(unsigned int) n)
		fatal_error("Can't read radius data",2);
}

void lj_sim::output_gnuplot(FILE *fp) {
	for(int i=0;i<n;i++) fprintf(fp,"%d %g %g %g %d\n",i,x[2*i],x[2*i+1],rad[i],ty[i]);
}

void lj_sim::output_povray(FILE *fp) {
	for(int i=0;i<n;i++) fprintf(fp,"sphere{<%.5f,%.5f,0>,%.5f}\n",x[2*i],x[2*i+1],rad[i]);
}
