#include "bi_interp.hh"
#include "shear_sim.hh"

bicubic_interp::bicubic_interp(shear_sim &ss_,int mode_) : ss(ss_),
    cor((1<<mode_)&(3|384)), mode(mode_), m(ss.gm), n(ss.gn),
    ax(ss.ax-(cor?2:1.5)*ss.dx), ay(ss.ay-(cor?2:1.5)*ss.dy), xsp(ss.xsp),
    ysp(ss.ysp), ijc(-1) {
    *a=a[1]=a[2]=a[3]=0;
    a[4]=a[5]=a[6]=a[7]=0;
    a[8]=a[9]=a[10]=a[11]=0;
    a[12]=a[13]=a[14]=a[15]=0;
}

void bicubic_interp::grid_index(double &x,double &y) {
    x=(x-ax)*xsp;y=(y-ay)*ysp;
    int i=static_cast<int>(x),j=static_cast<int>(y),ij;
    if(i<1||i>m-2||j<1||j>n-2) {
        fprintf(stderr,"Grid square (%d,%d) out of range in bicubic interpolation [%g,%g,%d]\n",i,j,x,y,mode);
        exit(1);
    }
    ij=i+m*j;
    if(ijc!=ij) table_setup(i,j,ij);
    x-=i;y-=j;
}

/** Sets up the table of coefficients of the bicubic interpolation function. */
void bicubic_interp::table_setup(int i,int j,int ij) {
    ijc=ij;
    c_field *up=ss.fbase+ij;
    double c00,c01,c02,c03;
    double c10,c11,c12,c13;
    double c20,c21,c22,c23;
    double c30,c31,c32,c33;
    compute_x(i,up,c01,c11,c21,c31);
    compute_x(i,up+m,c02,c12,c22,c32);
    compute_x(i,up-m,c00,c10,c20,c30);
    compute_x(i,up+2*m,c03,c13,c23,c33);
    fill_a(a,c00,c01,c02,c03);
    fill_a(a+4,c10,c11,c12,c13);
    fill_a(a+8,c20,c21,c22,c23);
    fill_a(a+12,c30,c31,c32,c33);
}

void bicubic_interp::compute_x(int i,c_field *up,double &c0,double &c1,double &c2,double &c3) {
    double f0=up[-1].fval(mode),f1=up->fval(mode),f2=up[1].fval(mode),f3=up[2].fval(mode);
    c0=f1;
    c1=-0.5*f0+0.5*f2;
    c2=f0-2.5*f1+2*f2-0.5*f3;
    c3=-0.5*f0+1.5*f1-1.5*f2+0.5*f3;
}

double bicubic_interp::f(double x,double y) {
    grid_index(x,y);
    return yl(a,y)+x*(yl(a+4,y)+x*(yl(a+8,y)+x*yl(a+12,y)));
}

double bicubic_interp::f_grad_f(double x,double y,double &fx,double &fy) {
    grid_index(x,y);
    fx=xsp*(yl(a+4,y)+x*(2*yl(a+8,y)+3*x*yl(a+12,y)));
    fy=ysp*(dyl(a,y)+x*(dyl(a+4,y)+x*(dyl(a+8,y)+x*dyl(a+12,y))));
    return yl(a,y)+x*(yl(a+4,y)+x*(yl(a+8,y)+x*yl(a+12,y)));
}

void bicubic_interp::fill_a(double *ap,double c0,double c1,double c2,double c3) {
    *ap=c1;
    ap[1]=-0.5*c0+0.5*c2;
    ap[2]=c0-2.5*c1+2*c2-0.5*c3;
    ap[3]=-0.5*c0+1.5*c1-1.5*c2+0.5*c3;
}
