#include <cmath>

#include "bd_sim.hh"
#include "extra_force.hh"
#include "qs_multi_bd.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent shear_sim
 * class.
 * \param[in] ss a reference to the parent shear_sim class. */
qs_multi_bd::qs_multi_bd(bd_sim &bd,int dom_li,int dom_ui,int dom_lj,int dom_uj) :
	m(dom_ui-dom_li+1), ml(bd.m), me(bd.me), n(dom_uj-dom_lj+1), mn(m*n),
	x_prd(false), y_prd(false), fm(bd.fm+dom_li+dom_lj*me), z(bd.vel),
	b(bd.src), eforce(bd.eforce), phi(bd.phi+dom_li+dom_lj*ml),
	c(bd.c+dom_li+dom_lj*ml), cc(bd.cc+dom_li+dom_lj*me), acc(1e-20),
	K(bd.K+bd.mu/3.0), qK(.25*K*bd.xsp*bd.ysp), Kmu(K+bd.mu), mu(bd.mu),
	dx(bd.dx), dy(bd.dy), xsp(bd.xsp), ysp(bd.ysp), xsp2(bd.xsp*bd.xsp),
	ysp2(bd.ysp*bd.ysp), xyfac(xsp2+ysp2), tmufac1(Kmu*xsp2+mu*ysp2),
	tmufac2(mu*xsp2+Kmu*ysp2), bcs(1), ocs(1e-6), ncs(0.14), rz(mat(0)),
	rc(mat(0,qK,qK,0)), mrc(mat(0,-qK,-qK,0)),
	d(new int[mn]), a_mem((m+n)*36), a(new mat[a_mem]) {}

qs_multi_bd::~qs_multi_bd() {
	delete [] a;
	delete [] d;
}

void qs_multi_bd::init(double visc,double dt) {
	na=0;

	// Set elasticity-related constants
	idt=1/dt;
	double viscdt=visc*idt;
	Kmx=(Kmu+viscdt)*xsp2;
	Kmy=(Kmu+viscdt)*ysp2;
	mx=(mu+viscdt)*xsp2;
	my=(mu+viscdt)*ysp2;
	cxx=-2*(tmufac1+viscdt*xyfac);
	cyy=-2*(tmufac2+viscdt*xyfac);

	// Set 2x2 matrices that are used within the multigrid solve
	rh=mat(Kmx,0,0,my);
	rv=mat(mx,0,0,Kmy);
	ru=mat(cxx,0,0,cyy);
	iu=mat(1/cxx,0,0,1/cyy);
	rb=mat(bcs*cxx,0,0,bcs*cyy);
	ib=mat(1/(bcs*cxx),0,0,1/(bcs*cyy));
	ro=mat(ocs*cxx,0,0,ocs*cyy);
	io=mat(1/(ocs*cxx),0,0,1/(ocs*cyy));
	rn=mat(ncs*cxx,0,0,ncs*cyy);
	in=mat(1/(ncs*cxx),0,0,1/(ncs*cyy));

	// Set up direction table
	double xfac=0.5/dt*xsp,yfac=0.5/dt*ysp,xbc=bcs*cxx,ybc=bcs*cyy;
#pragma omp parallel for
	for(int j=0;j<n;j++) {
	    double ul=fm[me*j].u,vl=fm[me*j].v;
	    double ur=fm[me*j+m-1].u,vr=fm[me*j+m-1].v;
	    for(int i=0,ij=j*m;i<m;i++,ij++) {

		// Deal with edges
		if(i==0||i==m-1) {
			c_field &f=fm[i+me*j];
			d[ij]=-2;
			b[ij].x=xbc*f.u;
			b[ij].y=ybc*f.v;
			continue;
		}

		// Deal with interior
		int *ccp=cc+i+me*j;
		if(*ccp<4) {
			d[ij]=-1;
			c_field *fp=fm+i+me*j,&f0=fp[-me-1],&f1=fp[-me],&f2=fp[-1],&f3=*fp;
			b[ij].x=(-f0.p+f1.p-f2.p+f3.p-f0.q+f1.q-f2.q+f3.q+f0.s-f1.s+f2.s-f3.s)*xfac
			       +(f0.tau+f1.tau-f2.tau-f3.tau)*yfac;
			b[ij].y=(f0.tau-f1.tau+f2.tau-f3.tau)*xfac
			       +(-f0.p-f1.p+f2.p+f3.p-f0.q-f1.q+f2.q+f3.q-f0.s-f1.s+f2.s+f3.s)*yfac;
			if(eforce!=NULL) b[ij]+=idt*eforce->force(i,j);
			continue;
		}

		// Deal with the Neumann-type free boundary condition
		int ijl=i+ml*j,*cp=c+ijl;
		if(*cp<4||cp[-1]<4||cp[-ml]<4||cp[-me]<4||ccp[-me-1]<4||ccp[-me]<4
		  ||ccp[-me+1]<4||ccp[-1]<4||ccp[1]<4||ccp[me-1]<4||ccp[me]<4||ccp[me+1]<4) {

			// Compute the normal vector
			double *phip=phi+ijl,
				nx=xsp*(*phip-phip[-1]+phip[-ml]-phip[-me]),
				ny=ysp*(*phip+phip[-1]-phip[-ml]-phip[-me]),
				nsq=nx*nx+ny*ny;
			if(nsq>1e-50) {
				nsq=1/sqrt(nsq);
				nx*=nsq;ny*=nsq;
			}

			// Store the coded direction indicator
			double fx=fabs(nx),fy=fabs(ny),z,zz;
			vec vsr;
			int k;
			mat h=xsp*mat(-Kmu*nx,-mu*ny,(mu-K)*ny,-mu*nx),
			    v=ysp*mat(-mu*ny,(mu-K)*nx,-mu*nx,-Kmu*ny),
			    dl(0),dc(0),dr(0),cl(0),cc(0),cr(0),ul(0),uc(0),ur(0);
			if(fx>fy) {
				zz=z=1-ny/fx;
				if(nx<0) {
					k=i;
					if(z<0.5) square1(z+0.5,h,v,dc,dr,cc,cr,true);
					else if (z<1.5) {z-=0.5;cc=-h;cr=h+v*(1-2*z);dr=-(1-z)*v;ur=z*v;}
					else square2(z-1.5,h,v,cc,cr,uc,ur,true);
				} else {
					k=i-1;
					if(z<0.5) square2(z+0.5,h,v,dl,dc,cl,cc,false);
					else if (z<1.5) {z-=0.5;cc=h;cl=-h+v*(1-2*z);dl=-(1-z)*v;ul=z*v;}
					else square1(z-1.5,h,v,cl,cc,ul,uc,false);
				}
				zz=f(zz);
				vsr=0.5*((2-zz)*contrib(k,j-1,nx,ny)+zz*contrib(k,j,nx,ny));
			} else {
				zz=z=1-nx/fy;
				if(ny<0) {
					k=j;
					if(z<0.5) square1(z+0.5,h,v,cl,cc,ul,uc,false);
					else if (z<1.5) {z-=0.5;cc=-v;uc=v+h*(1-2*z);ul=-(1-z)*h;ur=z*h;}
					else square2(2.5-z,h,v,cc,cr,uc,ur,true);
				} else {
					k=j-1;
					if(z<0.5) square2(0.5-z,h,v,dl,dc,cl,cc,false);
					else if (z<1.5) {z-=0.5;cc=v;dc=-v+h*(1-2*z);dl=-(1-z)*h;dr=z*h;}
					else square1(z-1.5,h,v,dc,dr,cc,cr,true);
				}
				zz=f(zz);
				vsr=0.5*((2-zz)*contrib(i-1,k,nx,ny)+zz*contrib(i,k,nx,ny));
			}

			// Store the source term
			cc=ru*(ncs/cc);
			b[ij]=idt*cc*vsr;

			// Precondition the elements
			dl=cc*dl;dc=cc*dc;dr=cc*dr;
			cl=cc*cl;cr=cc*cr;
			ul=cc*ul;uc=cc*uc;ur=cc*ur;

			// Store the matrix elements
#pragma omp critical
			{
				if(na==a_mem) add_normals_memory();
				d[ij]=na;
				a[na++]=dl;a[na++]=dc;a[na++]=dr;
				a[na++]=cl;a[na++]=cr;
				a[na++]=ul;a[na++]=uc;a[na++]=ur;
			}
		} else {
			d[ij]=-3;
			b[ij].x=xbc*(ul+double(i)/double(m-1)*(ur-ul));
			b[ij].y=ybc*(vl+double(i)/double(m-1)*(vr-vl));
		//	b[ij]=vec(0);
		}
	    }
	}

//	for(int j=0;j<n;j++) for(int i=0,ij=j*m;i<m;i++,ij++)
//		printf("%d %d %d dlook\n",i,j,d[ij]);
}

void qs_multi_bd::square1(double z,mat &h,mat &v,mat &a0,mat &a1,mat &a2,mat &a3,bool pos) {
	double fz=1-z;
	a0=-fz*(h+v);
	a1=fz*h-z*v;
	a2=fz*v-z*h;
	a3=z*(h+v);

	mat co=8*z*fz*(pos?a1:a2);
	a1-=co;a2-=co;a0+=co;a3+=co;
}

void qs_multi_bd::square2(double z,mat &h,mat &v,mat &a0,mat &a1,mat &a2,mat &a3,bool pos) {
	double fz=1-z;
	a0=-fz*h-z*v;
	a1=fz*(h-v);
	a2=z*(v-h);
	a3=fz*v+z*h;

	mat co=8*z*fz*(pos?a3:a0);
	a0-=co;a3-=co;a1+=co;a2+=co;
}

vec qs_multi_bd::contrib(int i,int j,double nx,double ny) {
	c_field &f=fm[i+me*j];
	double t0=f.p+f.q;
	return vec((-t0+f.s)*nx+f.tau*ny,f.tau*nx-(t0+f.s)*ny);
}

void qs_multi_bd::add_normals_memory() {
	printf("Double %d->%d\n",a_mem,a_mem<<1);
	mat *aa=new mat[a_mem<<1],*np=aa,*ne=a+na,*no=a;
	while(np<ne) *(np++)=*(no++);
	delete [] a;
	a=aa;
	a_mem<<=1;
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<qs_multi_bd,vec,mat>;
template class tgmg_base<tgmg_level<vec,mat>,vec,mat>;
template void tgmg_base<qs_multi_bd,vec,mat>::gauss_seidel();
template void tgmg_base<qs_multi_bd,vec,mat>::clear_z();
template double tgmg_base<qs_multi_bd,vec,mat>::mds();
template void tgmg_base<qs_multi_bd,vec,mat>::output(const char*,vec*,double,double,double,double);
