#include <cstdlib>

#include "grating.hh"

h_grating::h_grating(const int m_,const int n_,const double ax_,const double bx_,
			const double ay_,const double by_,const double k)
		: m(m_), n(n_), mn(m_*n_), x_prd(true), y_prd(false),
		acc(1e-20), ax(ax_), bx(bx_), ay(ay_), by(by_),
		dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)),
		kf(new double[mn-m]), A(new cpx[10*mn]),
		b(new cpx[11*mn]), x(b+mn), mg(*this,b,x) {}

/** The class destructor frees the dynamically allocated memory. */
h_grating::~h_grating() {
	delete [] b;
	delete [] A;
	delete [] kf;
}

/** Solves the Helmholtz problem using the preconditioned Bi-CGSTAB
 * method. */
void h_grating::solve() {
	int ij,l=0;
	cpx rho=1,alpha=1,omega=1,nrho,beta;
	cpx *p=b+2*mn,*v=b+3*mn,*rhat=b+4*mn,*r=b+5*mn,*s=b+6*mn,
	    *zz=b+7*mn,*kt=b+8*mn,*y=b+9*mn,*t=b+10*mn;

	// Set up the multigrid hierarchy
	mg.setup();

	// Set up
	for(ij=0;ij<mn;ij++) v[ij]=p[ij]=0;

	// Choose the arbitrary vector rhat to just be equal to the initial r
	mul_fa(x,r);
	for(ij=0;ij<mn;ij++) rhat[ij]=r[ij]=b[ij]-r[ij];

	// Do the Bi-CGSTAB iteration
	printf("Iter %d, L2 error %.10g\n",l,l2_error(r));
	while(true) {
		nrho=iprod(rhat,r);
		beta=(nrho/rho)*(alpha/omega);
		rho=nrho;
		for(ij=0;ij<mn;ij++) p[ij]=r[ij]+beta*(p[ij]-omega*v[ij]);
		pc_solve(p,y);
		mul_fa(y,v);
		alpha=rho/iprod(rhat,v);
		for(ij=0;ij<mn;ij++) s[ij]=r[ij]-alpha*v[ij];
		pc_solve(s,zz);
		mul_fa(zz,t);
		pc_solve(t,kt);
		omega=iprod(kt,zz)/iprod(kt,kt);
		for(ij=0;ij<mn;ij++) x[ij]=x[ij]+alpha*y[ij]+omega*zz[ij];
		l++;
		if(l==400) return;
		for(ij=0;ij<mn;ij++) r[ij]=s[ij]-omega*t[ij];
		printf("Iter %d, L2 error %.10g\n",l,l2_error(r));
	}
}

/** Calculates the inner product of two complex fields.
 * \param[in] (u,v) pointers to the two fields. */
cpx h_grating::iprod(cpx *u,cpx *v) {
	cpx w=std::conj(*u)*(*v);
	for(int ij=1;ij<mn;ij++) w+=std::conj(u[ij])*v[ij];
	return w;
}

/** Solves the preconditioning problem using multigrid.
 * \param[in] u a pointer to the source data.
 * \param[in] v a pointer to the solution. */
void h_grating::pc_solve(cpx *u,cpx *v) {
	mg.z=z=v;
	mg.mg[0]->y=v;
	mg.b=u;
	mg.solve_v_cycle();
}

/** Muultiplies an input field by the full problem A.
 * \param[in] u a pointer to the input field to consider.
 * \param[in] v a pointer to the solution. */
void h_grating::mul_fa(cpx *u,cpx *v) {
	int i,j,ij;
	z=u;
	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++)
		v[ij]=A[10*ij+4]*z[ij]+mul_a(i,ij);
}

/** Calculates the mean squared error of a complex field.
 * \param[in] u a pointer to the field.
 * \return The mean squared error. */
double h_grating::l2_error(cpx *u) {
	double l2=std::norm(*u);
	for(int ij=1;ij<mn;ij++) l2+=std::norm(u[ij]);
	return l2/mn;
}

void h_grating::setup_linear_system() {
	const double dydx=dy/dx,dxdy=dx/dy,fm(1./3.*(dxdy+dydx)),
	      fey(1./6.*(-2*dxdy+dydx)),fex(1./6.*(-2*dydx+dxdy)),
	      fc(-1./6.*(dxdy+dydx)),dxy=dx*dy,pc=dxy/36.,pex=dxy/18.,
	      pey=dxy/18.,pm=dxy/9.;
	const double alpha=1.;

	int i,j,ij;
	cpx *Ap;
	double *kfp,*kfl,k1,k2;

	for(ij=j=0;j<n;j++) {
		for(i=0;i<m;i++,ij++) {
			printf("%d %d %d\n",i,j,ij);
			Ap=A+10*ij;
			kfp=kf+ij;
			kfl=i==0?kfp+(m-1):kfp-1;

			for(int k=0;k<9;k++) Ap[k]=0;

			if(j==0) {
				Ap[4]-=4*fm;
				Ap[9]=1./Ap[4];
				continue;
			}

			if(j>0) {

				// Half of Laplacian stencil
				*Ap-=fc;Ap[1]-=2*fey;Ap[2]-=fc;
				Ap[3]-=fex;Ap[4]-=2*fm;Ap[5]-=fex;

				// Half of central stencil
				k1=kfl[-m];k2=kfp[-m];
				*Ap+=k1*pc;Ap[1]+=(k1+k2)*pey;Ap[2]+=k2*pc;
				Ap[3]+=k1*pex;Ap[4]+=(k1+k2)*pm;Ap[5]+=k2*pex;
			}

			if(j<n-1) {

				// Half of Laplacian stencil
				Ap[3]-=fex;Ap[4]-=2*fm;Ap[5]-=fex;
				Ap[6]-=fc;Ap[7]-=2*fey;Ap[8]-=fc;

				// Half of central stencil
				k1=*kfl;k2=*kfp;
				Ap[3]+=k1*pex;Ap[4]+=(k1+k2)*pm;Ap[5]+=k2*pex;
				Ap[6]+=k1*pc;Ap[7]+=(k1+k2)*pey;Ap[8]+=k2*pc;
			}

			if(j==n-1) {

				Ap[3]+=dx*alpha/6.;
				Ap[4]+=dx*alpha*(2/3.);
				Ap[5]+=dx*alpha/6.;
				b[ij]=dx;
			}
		//	for(int k=0;k<9;k++) printf("%d %g\n",k,Ap[k]);
			Ap[9]=1./Ap[4];
		}
	}
}

void h_grating::setup_test_problem() {
	int i,j,ij;
	double xx,yy;
	for(ij=j=0,yy=ay;j<n;j++,yy+=dy) {
		for(i=0,xx=ax;i<m;i++,xx+=dx,ij++) {
			x[ij]=0.;
			kf[ij]=-1;//j>32&&j<96?...:...;
			b[ij]=0.;
		}
	}
}

void h_grating::output(const char *filename,cpx *ff,int mode) {

	// Assemble the output filename and open the output file
	float *buf=new float[m+1];
	FILE *outf=fopen(filename,"wb");
	if(outf==NULL) {
		fputs("File output error\n",stderr);
		exit(1);
	}

	// Output the first line of the file
	int i,j;
	*buf=m;
	float *bp=buf+1,*be=buf+m+1;
	for(i=0;i<m;i++) *(bp++)=ax+i*dx;
	fwrite(buf,sizeof(float),m+1,outf);

	// Output the field values to the file
	cpx *fp=ff;
	for(j=0;j<n;j++) {
		*buf=ay+j*dy;bp=buf+1;
		switch(mode) {
			case 0:while(bp<be) *(bp++)=std::abs(*(fp++));break;
			case 1:while(bp<be) *(bp++)=std::arg(*(fp++));break;
			case 2:while(bp<be) *(bp++)=std::real(*(fp++));break;
			case 3:while(bp<be) *(bp++)=std::imag(*(fp++));
		}
		fwrite(buf,sizeof(float),m+1,outf);
	}

	// Close the file
	fclose(outf);
	delete [] buf;
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<h_grating,cpx,cpx>;
