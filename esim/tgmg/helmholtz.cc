#include <cstdlib>

const double pi=3.1415926535897932384626433832795;

#include "helmholtz.hh"

helmholtz::helmholtz(const int m_,const int n_,const double ax_,const double bx_,
			const double ay_,const double by_,const double k)
		: m(m_), n(n_), mn(m_*n_), x_prd(true), y_prd(false),
		acc(1e-20), ax(ax_), bx(bx_), ay(ay_), by(by_),
		dx((bx-ax)/(m-1)), dy((by-ay)/(n-1)),
		fm(cpx(-2/(dx*dx)-2/(dy*dy),k*k)), fm_inv(1./fm),
		fm_full(-2/(dx*dx)-2/(dy*dy)+k*k), fex(1./(dx*dx)),
		fey(1./(dy*dy)), b(new cpx[11*mn]), x(b+mn),
		mg(*this,b,x) {}

/** The class destructor frees the dynamically allocated memory. */
helmholtz::~helmholtz() {
	delete [] b;
}

/** Solves the Helmholtz problem using the preconditioned Bi-CGSTAB
 * method. */
void helmholtz::solve() {
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
cpx helmholtz::iprod(cpx *u,cpx *v) {
	cpx w=std::conj(*u)*(*v);
	for(int ij=1;ij<mn;ij++) w+=std::conj(u[ij])*v[ij];
	return w;
}

/** Solves the preconditioning problem using multigrid.
 * \param[in] u a pointer to the source data.
 * \param[in] v a pointer to the solution. */
void helmholtz::pc_solve(cpx *u,cpx *v) {
	mg.z=z=v;
	mg.mg[0]->y=v;
	mg.b=u;
	mg.solve_v_cycle();
}

/** Muultiplies an input field by the full problem A.
 * \param[in] u a pointer to the input field to consider.
 * \param[in] v a pointer to the solution. */
void helmholtz::mul_fa(cpx *u,cpx *v) {
	int i,j,ij;
	z=u;
	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++)
		v[ij]=fm_full*z[ij]+mul_a(i,ij);
}

/** Calculates the mean squared error of a complex field.
 * \param[in] u a pointer to the field.
 * \return The mean squared error. */
double helmholtz::l2_error(cpx *u) {
	double l2=std::norm(*u);
	for(int ij=1;ij<mn;ij++) l2+=std::norm(u[ij]);
	return l2/mn;
}

void helmholtz::setup_test_problem() {
	int i,j,ij;
	double xx,yy;
	for(ij=j=0,yy=ay;j<n;j++,yy+=dy) {
		for(i=0,xx=ax;i<m;i++,xx+=dx,ij++) {
			x[ij]=0.;
			b[ij]=0.;
		}
	}
	double fac=2*pi/(bx-ax);
	for(ij=mn-m,i=0,xx=ax;i<m;i++,xx+=dx,ij++) {
		b[ij]=sin(4*fac*xx)*fm_full;
	}
}

void helmholtz::output(const char *filename,cpx *ff,int mode) {

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
template class tgmg<helmholtz,cpx,cpx>;
template class tgmg_base<tgmg_level<cpx,cpx>,cpx,cpx>;
