#include "sbfrac.hh"

#define VISCOSITY

sbfrac::sbfrac(int im,int in,double iax,double ibx,double iay,double iby,
	     double ichi_inf,double ic_0_inv,double inu,double imu,double iK,double itmult,const char *ifilename)
	: m(im),n(in),mn(im*in),ax(iax),ay(iay),dx((ibx-iax)/(im-1)),dy((iby-iay)/(in-1)),
	xsp(1/dx),ysp(1/dy),chi_inf(ichi_inf),c_0_inv(ic_0_inv),nu(inu),mu(imu),K(iK),mu_inv(1/imu),tmult(itmult),filename(ifilename),
	ls(im,in,iax,iay,dx,dy,-5*sqrt(dx*dx+dy*dy),2.5*sqrt(dx*dx+dy*dy)),
	phi(ls.phi), c(ls.s) {
	u=new double[mn];v=new double[mn];p=new double[mn];
	s=new double[mn],tau=new double[mn];chi=new double[mn];
	cu=new double[mn];cv=new double[mn];cp=new double[mn];
	cs=new double[mn];ctau=new double[mn];cchi=new double[mn];

#ifdef TRACERS
	trace=new double[2*TRACERS];
#endif
	init_fields();

	ls.add_field(u);ls.add_field(v);ls.add_field(p);
	ls.add_field(s);ls.add_field(tau);ls.add_field(chi);
	ls.build_band();
}

sbfrac::~sbfrac() {
	delete [] u;delete [] v;delete [] p;
	delete [] s;delete [] tau;delete [] chi;
	delete [] cu;delete [] cv;delete [] cp;
	delete [] cs;delete [] ctau;delete [] cchi;
}

void sbfrac::relax_fields() {
	const double rela=0.17;
	const double relaz=0;
	double uu,vu,pu,su,tauu,chiu;
	double ud,vd,pd,sd,taud,chid;
	double px,py,sx,sy,taux,tauy;
	double fx,fy,res;
	int i,j,ij,l;
	for(ij=0;ij<mn;ij++) cs[ij]=ctau[ij]=0;
	for(l=0;l<300;l++) {
		ls.extrapolate_fields();
		res=0;
		for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) if(c[ij]<4) {
			compute_bc_deriv(uu,vu,pu,su,tauu,chiu,i,j,i+1,j);
			compute_bc_deriv(ud,vd,pd,sd,taud,chid,i,j,i-1,j);
			px=(pu-pd);sx=(su-sd);taux=(tauu-taud);
			cs[ij]+=(su-2*s[ij]+sd)*rela*0.01;

			compute_bc_deriv(uu,vu,pu,su,tauu,chiu,i,j,i+1,j);
			compute_bc_deriv(ud,vd,pd,sd,taud,chid,i,j,i-1,j);
			py=(pu-pd);sy=(su-sd);tauy=(tauu-taud);
			cs[ij]+=(su-2*s[ij]+sd)*rela*0.01;
		//	ctau[ij]+=(tauu-2*tau[ij]+taud)*rela*0.0001;

			fx=(-px+sx+tauy);
			fy=(-py-sy-taux);
			res+=fx*fx+fy*fy;

			if(i>0&&c[ij-1]<4) {cs[ij-1]+=rela*fx;ctau[ij-1]-=relaz*fy;}
			if(i<m-1&&c[ij+1]<4) {cs[ij+1]-=rela*fx;ctau[ij+1]+=relaz*fy;}
		//	if(j>0&&c[ij-m]<4) {cs[ij-m]-=relaz*fy;ctau[ij-m]+=rela*fx;}
		//	if(j<n-1&&c[ij+m]<4) {cs[ij+m]+=relaz*fy;ctau[ij+m]-=rela*fx;}
		}
		if(l%10==0) write_files(l);
		for(ij=0;ij<mn;ij++) if(c[ij]<4) {
			s[ij]+=cs[ij];cs[ij]=0;
			tau[ij]+=ctau[ij];ctau[ij]=0;
		}
		cout << l << " " << res << endl;
	}
}

void sbfrac::solve(double t_start,double t_end,int frames) {
	bool ou,od;
	time=t_start;
	double time_interval=(t_end-t_start)/frames,target_time;
	double hx,hy,dt;
	double u1,v1,p1,s1,tau1,chi1;
	double u2,v2,p2,s2,tau2,chi2;
	double ux,vx,px,sx,taux;
	double uy,vy,py,sy,tauy;
	double uu,vu,pu,su,tauu,chiu,ud,vd,pd,sd,taud,chid,tomega,magn,magr;
	const double kappa=0.5;
	double sT,sC,rho,rhoc,pM,M,xi,Gamma,dplas;
#ifdef VISCOSITY
	double hxx,hyy,uxx,vxx;
#endif

	bool output;
	int i,j,ij,k;

//	const double timestep=0.00005;
	static char fname[256];
	sprintf(fname,"%s/kfile",filename);
	kfile.open(fname,fstream::out|fstream::trunc);

	write_files(0);

	for(k=1;k<=frames;k++) {
		target_time=t_start+time_interval*k;
		output=true;
		while (output) {
			dt=dx*dx*tmult;
			if(time+dt+1e-11 >= target_time) {
				dt=target_time-time;
				output=false;
			}

			hx=0.5*dt*xsp;hy=0.5*dt*ysp;
#ifdef VISCOSITY
			hxx=2*hx*xsp;hyy=2*hy*ysp;
#endif
			ls.extrapolate_fields();

			// Update level set
			ls.move(this,dt);

			// Update time and boundary pressure
			time+=dt;

			// Update the fields
			ij=0;
			for(j=0;j<n;j++) for(i=0;i<m;i++) {
				ij=i+m*j;
				if(c[ij]<4) {

					// Compute x derivatives, taking into account BCs
					ou=compute_bc_deriv(uu,vu,pu,su,tauu,chiu,i,j,i+1,j);
					od=compute_bc_deriv(ud,vd,pd,sd,taud,chid,i,j,i-1,j);
					tomega=hx*(vu-vd);

					// Calculate a one-sided derivative using ENO2
					// if the gridpoints are available
					if(u[ij]>0) {
						if(od&&i>1&&c[ij-2]<4) {
							u1=hx*eno2(uu,u[ij],ud,u[ij-2]);v1=hx*eno2(vu,v[ij],vd,v[ij-2]);
							p1=hx*eno2(pu,p[ij],pd,p[ij-2]);s1=hx*eno2(su,s[ij],sd,s[ij-2]);
							tau1=hx*eno2(tauu,tau[ij],taud,tau[ij-2]);chi1=hx*eno2(chiu,chi[ij],chid,chi[ij-2]);
						} else {
							u1=2*hx*(u[ij]-ud);v1=2*hx*(v[ij]-vd);
							p1=2*hx*(p[ij]-pd);s1=2*hx*(s[ij]-sd);
							tau1=2*hx*(tau[ij]-taud);chi1=2*hx*(chi[ij]-chid);
						}
					} else {
						if(ou&&i<m-2&&c[ij+2]<4) {
							u1=-hx*eno2(ud,u[ij],uu,u[ij+2]);v1=-hx*eno2(vd,v[ij],vu,v[ij+2]);
							p1=-hx*eno2(pd,p[ij],pu,p[ij+2]);s1=-hx*eno2(sd,s[ij],su,s[ij+2]);
							tau1=-hx*eno2(taud,tau[ij],tauu,tau[ij+2]);chi1=-hx*eno2(chid,chi[ij],chiu,chi[ij+2]);
						} else {
							u1=2*hx*(uu-u[ij]);v1=2*hx*(vu-v[ij]);
							p1=2*hx*(pu-p[ij]);s1=2*hx*(su-s[ij]);
							tau1=2*hx*(tauu-tau[ij]);chi1=2*hx*(chiu-chi[ij]);
						}
					}
					ux=hx*(uu-ud);vx=hx*(vu-vd);px=hx*(pu-pd);sx=hx*(su-sd);taux=hx*(tauu-taud);
#ifdef VISCOSITY
					uxx=hxx*(uu-2*u[ij]+ud);
					vxx=hxx*(vu-2*v[ij]+vd);
#endif

					// Compute y derivatives, taking into account BCs
					ou=compute_bc_deriv(uu,vu,pu,su,tauu,chiu,i,j,i,j+1);
					od=compute_bc_deriv(ud,vd,pd,sd,taud,chid,i,j,i,j-1);
					tomega+=hy*(ud-uu);

					// Calculate a one-sided derivative using ENO2
					// if the gridpoints are available
					if(v[ij]>0) {
						if(od&&j>1&&c[ij-2*m]<4) {
							u2=hy*eno2(uu,u[ij],ud,u[ij-2*m]);v2=hy*eno2(vu,v[ij],vd,v[ij-2*m]);
							p2=hy*eno2(pu,p[ij],pd,p[ij-2*m]);s2=hy*eno2(su,s[ij],sd,s[ij-2*m]);
							tau2=hy*eno2(tauu,tau[ij],taud,tau[ij-2*m]);chi2=hy*eno2(chiu,chi[ij],chid,chi[ij-2*m]);
						} else {
							u2=2*hy*(u[ij]-ud);v2=2*hy*(v[ij]-vd);
							p2=2*hy*(p[ij]-pd);s2=2*hy*(s[ij]-sd);
							tau2=2*hy*(tau[ij]-taud);chi2=2*hy*(chi[ij]-chid);
						}
					} else {
						if(ou&&j<n-2&&c[ij+2*m]<4) {
							u2=-hy*eno2(ud,u[ij],uu,u[ij+2*m]);v2=-hy*eno2(vd,v[ij],vu,v[ij+2*m]);
							p2=-hy*eno2(pd,p[ij],pu,p[ij+2*m]);s2=-hy*eno2(sd,s[ij],su,s[ij+2*m]);
							tau2=-hy*eno2(taud,tau[ij],tauu,tau[ij+2*m]);chi2=-hy*eno2(chid,chi[ij],chiu,chi[ij+2*m]);
						} else {
							u2=2*hy*(uu-u[ij]);v2=2*hy*(vu-v[ij]);
							p2=2*hy*(pu-p[ij]);s2=2*hy*(su-s[ij]);
							tau2=2*hy*(tauu-tau[ij]);chi2=2*hy*(chiu-chi[ij]);
						}
					}
					uy=hy*(uu-ud);vy=hy*(vu-vd);py=hy*(pu-pd);sy=hy*(su-sd);tauy=hy*(tauu-taud);
					cu[ij]=-u[ij]*u1-v[ij]*u2+(-px+sx+tauy)*mu_inv;
					cv[ij]=-u[ij]*v1-v[ij]*v2+(-py-sy+taux)*mu_inv;
					cp[ij]=-u[ij]*p1-v[ij]*p2-K*(ux+vy);
					magr=sqrt(magn=s[ij]*s[ij]+tau[ij]*tau[ij]);

					sT=tanh(TE/T*sinh(magr/mub));
					sC=fsC(magr);
					rho=frho(T);
					if(magr>1e-10) {
						pM=1+magr*sT+rho/(2*sC);
						M=(pM-sqrt(pM*pM-4*magr*sT))/(2*magr);
					} else {
						pM=1+1e-10*sT+rho/(2*sC);
						M=(pM-sqrt(pM*pM-4e-10*sT))/(2e-10);
					}
					if(1-M>1e-10) {
						xi=mub*asinh((T/TE)*atanh(M));
					} else {
						rhoc=rho/(2*sC);
						if((magr-1)*(magr-1)>rhoc) {
							if(magr>1) {
								xi=atanh(1/magr+rhoc/(magr-magr*magr));
							} else {
		//						if(1-sT<1e-10) {
		//							xi=0.5*(log(1+sT-rhoc/(1-magr*sT))-log(rhoc*sT/(1-magr*sT)));
		//						} else {
									xi=0.5*(log(1+sT-rhoc/(1-magr*sT))-log((1-sT)+rhoc*sT/(1-magr*sT)));
		//						}
							//	xi=0.5*(log(1+sT+rhoc/(magr-1))-log(rhoc)+log(1-magr));
							}
						} else {
							xi=0.5*(log(2)-log((magr-1)*0.5-rhoc*0.5+0.5*sqrt((magr-1)*(magr-1)+2*rhoc*(2+(magr-1))+rhoc*rhoc)));
						}
						xi=mub*asinh(T/TE*xi);
					}
					Gamma=(2*sC*(sT-M)*(magr-xi)+M*xi*rho)/(1-M*xi);

					dplas=0.5*dt*nu*nu_z*c_0_inv*exp(-1/chi[ij])*chi[ij];
				//	cchi[ij]=chi[ij]-u[ij]*chi1-v[ij]*chi2+dt*nu*c_0_inv*exp(-1/chi[ij])*chi[ij]*(Gamma*(chi_inf-chi[ij])/chi_inf+kappa*rho*(1-TZ*chi[ij]/T));
					cchi[ij]=(chi[ij]-u[ij]*chi1-v[ij]*chi2+dplas*(Gamma+kappa*rho))/(1+dplas*(Gamma/chi_inf+kappa*rho*TZ/T));
					if(isnan(cchi[ij])) cout << ij << "messed" << endl;

					dplas=1/(1+dt*mu*nu*sC*exp(-1/chi[ij])*(sT-M));

					cs[ij]=(s[ij]-u[ij]*s1-v[ij]*s2-tomega*tau[ij]+mu*(ux-vy))*dplas;
					ctau[ij]=(tau[ij]-u[ij]*tau1-v[ij]*tau2+tomega*s[ij]+mu*(uy+vx))*dplas;
#ifdef VISCOSITY
					cu[ij]+=visc*(uxx+hyy*(uu-2*u[ij]+ud));
					cv[ij]+=visc*(vxx+hyy*(vu-2*v[ij]+vd));
#endif
				}
			}

#ifdef TRACERS
			for(i=0;i<2*TRACERS;i+=2) {
				u1=dt*ls.bicubic(u,trace[i],trace[i+1]);
				v1=dt*ls.bicubic(v,trace[i],trace[i+1]);
				trace[i]+=u1;
				trace[i+1]+=v1;
			}
#endif

			//Apply the changes
			for(j=0;j<n;j++) for(i=0;i<m;i++) {
				ij=i+m*j;
				if(c[ij]<4) {
					u[ij]+=cu[ij];cu[ij]=0;
					v[ij]+=cv[ij];cv[ij]=0;
					p[ij]+=cp[ij];cp[ij]=0;
					s[ij]=cs[ij];cs[ij]=0;
					tau[ij]=ctau[ij];ctau[ij]=0;
					chi[ij]=cchi[ij];cchi[ij]=0;
				}
			}
			if(!nucleated) nucleate_void();

		//	cout << "Time: " << time << endl;
			print_maxima();
		}

//		cout << time << endl;

		cout << "# Output frame " << k << endl;

		//Clean up grid before output
		write_files(k);
	}
	kfile.close();
}

void sbfrac::print_maxima() {
	double minu,maxu,minp,maxp,minchi,maxchi;
	int i,j,ij=0,maxi=-1;
	while(c[ij]>=4) ij++;
	minu=maxu=u[ij];
	minp=maxp=p[ij];
	minchi=maxchi=chi[ij];
	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) {
		if(c[ij]<4) {
			if(u[ij]<minu) minu=u[ij];
			if(u[ij]>maxu) maxu=u[ij];
			if(p[ij]<minp) minp=p[ij];
			if(p[ij]>maxp) maxp=p[ij];
			if(chi[ij]<minchi) minchi=chi[ij];
			if(chi[ij]>maxchi) maxchi=chi[ij];
		} else if(i>maxi) maxi=i;
	}
	printf("Ext: %f %f %f %f %f %f %f %f\n",time,minu,maxu,minp,maxp,minchi,maxchi,sinten*(lambda+pgamma*time));
}

void sbfrac::nucleate_void() {
//	if(!nucleated) {
//	} else {
//		pressure_cutoff=-2.7;void_size=0.5*dx;
//	}
	bool any_voids=false;
	int i,j,ij;
	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) if(c[ij]<4&&p[ij]<pressure_cutoff) {
		cout << "void at " << i << " " << j << endl;
		int i2,j2,ij2;double xs,ys,phit;
		for(ij2=j2=0;j2<n;j2++) for(i2=0;i2<m;i2++,ij2++) {
			xs=dx*(i2-i);ys=dy*(j2-j);
			phit=void_size*dx-sqrt(xs*xs+ys*ys);
			if(phit>phi[ij2]) phi[ij2]=phit;
			any_voids=true;
		//	nucleated=true;
		}
	}
	if(any_voids) ls.build_band();
}

inline double sbfrac::eno2(double p0,double p1,double p2,double p3) {
	return abs(p0-2*p1+p2)>abs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

inline double sbfrac::qs(double s) {
	return fabs(s)>1?(fabs(s)-1)*(fabs(s)-1):0;
}

inline double sbfrac::qs_new(double s,double cchi) {
	s=fabs(s);
	double temp=tanh(s/cchi)-1/s;
	return temp>0?cosh(s/cchi)*temp:0;
}

bool sbfrac::compute_bc_deriv(double &uu,double &vu,double &pu,double &su,double &tauu,double &chiu,int i,int j,int ni,int nj) {
	if (ni<0||ni>=m||nj<0||nj>=n) {
		int ij=i+m*j,bij=2*i-ni+m*(2*j-nj);
		double x=ax+dx*ni,y=ay+dy*nj;
		double r=sqrt(x*x+y*y),theta=arg(x,y);
		uu=Factor*sqrt(r)*(mode_ii?(2*Kappa+3)*sin(theta*0.5)+sin(theta*1.5)
				:(2*Kappa-1)*cos(theta*0.5)-cos(theta*1.5))+drift;
		vu=Factor*sqrt(r)*(mode_ii?-(2*Kappa-3)*cos(theta*0.5)-cos(theta*1.5)
				:(2*Kappa+1)*sin(theta*0.5)-sin(theta*1.5));
		pu=2*p[ij]-p[bij];
		su=2*s[ij]-s[bij];
		tauu=2*tau[ij]-tau[bij];
		chiu=2*chi[ij]-chi[bij];
		return false;
	} else {
		int nij=ni+m*nj;
		uu=u[nij];vu=v[nij];chiu=chi[nij];
		if (c[nij]<4) {
			pu=p[nij];su=s[nij];tauu=tau[nij];
			return true;
		} else {
			int ij=i+m*j;
			double sigmap,magn,n1p,n2p,sin2t,cos2t,z=phi[ij]/(phi[ij]-phi[nij]),
			pp=p[ij]*(1-z)+z*p[nij],
			sp=s[ij]*(1-z)+z*s[nij],
			taup=tau[ij]*(1-z)+z*tau[nij];
			n1p=phi_x(i,ij)*(1-z)+z*phi_x(ni,nij);
			n2p=phi_y(i,ij)*(1-z)+z*phi_y(ni,nij);
			magn=sqrt(n1p*n1p+n2p*n2p);
			if(magn>1e-8) {
				n1p/=magn;n2p/=magn;
				sin2t=2*n1p*n2p;
				cos2t=1-2*n2p*n2p;
				sigmap=-pp-sp*cos2t-taup*sin2t;
				pp=-0.5*sigmap;
				sp=-0.5*sigmap*cos2t;
				taup=-0.5*sigmap*sin2t;
				int bi=2*i-ni,bj=2*j-nj;
				if(bi>=0&&bi<m&&bj>=0&&bj<n) {
					int bij=bi+m*bj;
					pu=(2*pp-p[bij]*(1-z))/(1+z);
					su=(2*sp-s[bij]*(1-z))/(1+z);
					tauu=(2*taup-tau[bij]*(1-z))/(1+z);
				} else {
					printf("#Outside elt %d %d %d %d %g\n",i,j,ni,nj,z);
					pu=(pp-p[ij]*(1-z))/z;
					su=(sp-s[ij]*(1-z))/z;
					tauu=(taup-tau[ij]*(1-z))/z;
				}
			} else {
				pu=p[nij];su=s[nij];tauu=tau[nij];
			}
			return false;
		}
	}
}

inline double sbfrac::phi_x(int i,int ij) {
	return (i==0?phi[ij+1]-phi[ij]:(i==m-1?phi[ij]-phi[ij-1]:(phi[ij+1]-phi[ij-1])*0.5))*xsp;
}

inline double sbfrac::phi_y(int i,int ij) {
	return (ij<m?phi[ij+m]-phi[ij]:(ij>=mn-m?phi[ij]-phi[ij-m]:(phi[ij+m]-phi[ij-m])*0.5))*ysp;
}

void sbfrac::clean() {
	for(int ij=0;ij<mn;ij++) if(c[ij]>3) u[ij]=v[ij]=p[ij]=s[ij]=tau[ij]=chi[ij]=0;
}

void sbfrac::write_files(int k) {
	kfile << k << " " << time << " " << sinten*(lambda+gamma*time) << endl;

	output(filename,"u",u,k);
	output(filename,"v",v,k);
	output(filename,"p",p,k);
	output(filename,"s",s,k);
	//output(filename,"tau",tau,k);
	output(filename,"chi",chi,k);
	//output(filename,"c",c,k);
	output(filename,"phi",phi,k);
	for(int ij=0;ij<mn;ij++) cu[ij]=sqrt(s[ij]*s[ij]+tau[ij]*tau[ij]);
	output(filename,"dev",cu,k);
#ifdef TRACERS
	ofstream outfile;
	static char fname[256];sprintf(fname,"%s/trace.%d",filename,k);
	outfile.open(fname,fstream::out|fstream::trunc);
	for(int i=0;i<2*TRACERS;i+=2) outfile << trace[i] << " " << trace[i+1] << "\n";
	outfile.close();
#endif
}

template<class T>
void sbfrac::output(const char* filename,const char *suffix,T *array,const int sn) {
	int i,j,ij=0,ds=sizeof(float);
	float data_float;const char *pfloat;
	pfloat=(const char*)&data_float;

	ofstream outfile;
	static char fname[256];sprintf(fname,"%s/%s.%d",filename,suffix,sn);
	outfile.open(fname,fstream::out|fstream::trunc|fstream::binary);

	data_float=m;outfile.write(pfloat,ds);
	for(i=0;i<m;i++) {
		data_float=ax+i*dx;
		outfile.write(pfloat,ds);
	}
	for(ij=j=0;j<n;j++) {
		data_float=ay+j*dy;
		outfile.write(pfloat,ds);
		for(i=0;i<m;i++,ij++) {
			data_float=array[ij];
			outfile.write(pfloat,ds);
		}
	}

	outfile.close();
}

double sbfrac::velocity(double x,double y) {
	double x0,y0;
	if(x<-0.99) x=-0.99;
	ls.bicubic(x,y,x0,y0);
//	double tm2=pressure_decay-bicubic(ppp,x,y);
//	tm2=tm2>0?tm2*12:0;
	double temp=-ls.bicubic(u,x,y)*x0-ls.bicubic(v,x,y)*y0;//+tm2;
	return temp;
}

void sbfrac::fatal_error(const char *p,int code) {
	cerr << p << endl;
	cerr << "Current time: " << time << endl;
	exit(code);
}

inline double sbfrac::arg(double x,double y) {
	return x+y>0?(x>y?atan(y/x):0.5*pi-atan(x/y)):(x>y?-atan(x/y)-0.5*pi:atan(y/x)+(y>0?pi:-pi));
}

inline double sbfrac::frho(double T) {
	const double a=3;
	return T>T0?exp(-T1/(T-T0)*exp(-a*(T-T0)/(TA-T0))):0;
}

inline double sbfrac::fsC(double s) {
	const double s1=1.0/1.1;
	double fac=(TE/T)*cosh(s/mub);
	return (fac>100?0.5*exp(-(TE/T)*exp(-s/mub)):exp(-fac)*cosh((TE/T)*sinh(s/mub)))*sqrt(1+(s/s1)*(s/s1));
}
