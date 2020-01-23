#include "sbsim.hh"
using namespace std;

sbsim::sbsim(int im,int in,double iax,double ibx,double iay,double iby,
	     double ichi_inf,double ic_0_inv,double inu,double imu,double iK,double itmult,const char *ifilename)
	: m(im),n(in),mn(im*in),ax(iax),ay(iay),dx((ibx-iax)/(im-1)),dy((iby-iay)/(in-1)),
	xsp(1/dx),ysp(1/dy),chi_inf(ichi_inf),c_0_inv(ic_0_inv),nu(inu),mu(imu),K(iK),mu_inv(1/imu),tmult(itmult),
	filename(ifilename),
	ls(im,in,iax,iay,dx,dy,-5*sqrt(dx*dx+dy*dy),2.5*sqrt(dx*dx+dy*dy)),
	phi(ls.phi), c(ls.s) {
	u=new double[mn];v=new double[mn];p=new double[mn];
	s=new double[mn],tau=new double[mn];chi=new double[mn];
	cu=new double[mn];cv=new double[mn];cp=new double[mn];
	cs=new double[mn];ctau=new double[mn];cchi=new double[mn];

	init_fields();

	e_reg.setup_fields(6,u,v,p,s,tau,chi);
	ls.build_band();
}

sbsim::~sbsim() {
	delete [] u;delete [] v;delete [] p;
	delete [] s;delete [] tau;delete [] chi;
	delete [] cu;delete [] cv;delete [] cp;
	delete [] cs;delete [] ctau;delete [] cchi;
}

void sbsim::solve(double t_start,double t_end,int frames) {
	bool ou,od;
	time=t_start;
	double time_interval=(t_end-t_start)/frames,target_time;
	double hx,hy,dt;
	double u1,v1,p1,s1,tau1,chi1;
	double u2,v2,p2,s2,tau2,chi2;
#ifdef FULL_ENO
	const double Kmf=sqrt(1+K*mu_inv);
	double chix,chiy,t1,t2,t3,t4,t5;
#endif
	double ux,vx,px,sx,taux;
	double uy,vy,py,sy,tauy;
	double uu,vu,pu,su,tauu,chiu,ud,vd,pd,sd,taud,chid,tomega,magn,magr,cha;
#ifdef VISCOSITY
	const double visc=VISCOSITY;
	double hxx,hyy,uxx,vxx;
#endif

	min_i=int((-wallx-ax)*xsp)+1;if(min_i<0) min_i=0;
	max_i=int((wallx-ax)*xsp);if(max_i>=m) max_i=m-1;

	bool output;
	int i,j,ij,k,l;

//	const double timestep=0.00005;
	static char fname[256];
	sprintf(fname,"%s/wall",filename);
	wallfile.open(fname,fstream::out|fstream::trunc|fstream::binary);
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
			simple_extrapolation(min_i-1,min_i,min_i+1);
			simple_extrapolation(max_i+1,max_i,max_i-1);
			ls.extrapolate_fields(e_reg);

			// Update level set
			ls.move(this,dt);

			// Update time and boundary pressure
			time+=dt;

			wallx+=wallu*dt;
			l=int((-wallx-ax)*xsp)+1;if(l<0) l=0;
			if(l<min_i) simple_extrapolation(l,l+1,l+2);
			min_i=l;

			l=int((wallx-ax)*xsp);if(l>=m) l=m-1;
			if(l>max_i) simple_extrapolation(l-2,l-1,l);
			max_i=l;

			// Update the fields
			ij=0;
			for(j=0;j<n;j++) for(i=min_i;i<=max_i;i++) {
				ij=i+m*j;
				if(c[ij]<4) {

					// Compute x derivatives, taking into account BCs
					ou=compute_bc_deriv(uu,vu,pu,su,tauu,chiu,i,j,i+1,j);
					od=compute_bc_deriv(ud,vd,pd,sd,taud,chid,i,j,i-1,j);
					tomega=hx*(vu-vd);

					// Calculate a one-sided derivative using ENO2
					// if the gridpoints are available
#ifdef FULL_ENO
					if(od&&i>1&&c[ij-2]<4) {
						u1=hx*eno2(uu,u[ij],ud,u[ij-2]);v1=hx*eno2(vu,v[ij],vd,v[ij-2]);
						p1=hx*eno2(pu,p[ij],pd,p[ij-2]);s1=hx*eno2(su,s[ij],sd,s[ij-2]);
						tau1=hx*eno2(tauu,tau[ij],taud,tau[ij-2]);chi1=hx*eno2(chiu,chi[ij],chid,chi[ij-2]);
					} else {
						u1=2*hx*(u[ij]-ud);v1=2*hx*(v[ij]-vd);
						p1=2*hx*(p[ij]-pd);s1=2*hx*(s[ij]-sd);
						tau1=2*hx*(tau[ij]-taud);chi1=2*hx*(chi[ij]-chid);
					}
					if(ou&&i<m-2&&c[ij+2]<4) {
						u2=-hx*eno2(ud,u[ij],uu,u[ij+2]);v2=-hx*eno2(vd,v[ij],vu,v[ij+2]);
						p2=-hx*eno2(pd,p[ij],pu,p[ij+2]);s2=-hx*eno2(sd,s[ij],su,s[ij+2]);
						tau2=-hx*eno2(taud,tau[ij],tauu,tau[ij+2]);chi2=-hx*eno2(chid,chi[ij],chiu,chi[ij+2]);
					} else {
						u2=2*hx*(uu-u[ij]);v2=2*hx*(vu-v[ij]);
						p2=2*hx*(pu-p[ij]);s2=2*hx*(su-s[ij]);
						tau2=2*hx*(tauu-tau[ij]);chi2=2*hx*(chiu-chi[ij]);
					}
					chix=u[ij]>0?chi1:chi2;
					t1=u[ij]>0?mu*p1+s1*K:mu*p2+s2*K;
					t2=u[ij]-1>0?mu*v1+tau1:mu*v2+tau2;
					t3=u[ij]+1>0?-mu*v1+tau1:-mu*v2+tau2;
					t4=u[ij]-Kmf>0?mu*u1*Kmf-p1+s1:mu*u2*Kmf-p2+s2;
					t5=u[ij]+Kmf>0?-mu*u1*Kmf-p1+s1:-mu*u2*Kmf-p2+s2;

					ux=0.5*mu_inv/Kmf*(t4-t5);
					vx=0.5*mu_inv*(t2-t3);
					px=(t1-0.5*K*(t4+t5))/(K+mu);
					sx=(t1+0.5*mu*(t4+t5))/(K+mu);
					taux=0.5*(t2+t3);
#else
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
#endif
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
#ifdef FULL_ENO
					if(od&&j>1&&c[ij-2*m]<4) {
						u1=hy*eno2(uu,u[ij],ud,u[ij-2*m]);v1=hy*eno2(vu,v[ij],vd,v[ij-2*m]);
						p1=hy*eno2(pu,p[ij],pd,p[ij-2*m]);s1=hy*eno2(su,s[ij],sd,s[ij-2*m]);
						tau1=hy*eno2(tauu,tau[ij],taud,tau[ij-2*m]);chi1=hy*eno2(chiu,chi[ij],chid,chi[ij-2*m]);
					} else {
						u1=2*hy*(u[ij]-ud);v1=2*hy*(v[ij]-vd);
						p1=2*hy*(p[ij]-pd);s1=2*hy*(s[ij]-sd);
						tau1=2*hy*(tau[ij]-taud);chi1=2*hy*(chi[ij]-chid);
					}
					if(ou&&j<n-2&&c[ij+2*m]<4) {
						u2=-hy*eno2(ud,u[ij],uu,u[ij+2*m]);v2=-hy*eno2(vd,v[ij],vu,v[ij+2*m]);
						p2=-hy*eno2(pd,p[ij],pu,p[ij+2*m]);s2=-hy*eno2(sd,s[ij],su,s[ij+2*m]);
						tau2=-hy*eno2(taud,tau[ij],tauu,tau[ij+2*m]);chi2=-hy*eno2(chid,chi[ij],chiu,chi[ij+2*m]);
					} else {
						u2=2*hy*(uu-u[ij]);v2=2*hy*(vu-v[ij]);
						p2=2*hy*(pu-p[ij]);s2=2*hy*(su-s[ij]);
						tau2=2*hy*(tauu-tau[ij]);chi2=2*hy*(chiu-chi[ij]);
					}
					chiy=v[ij]>0?chi1:chi2;
					t1=v[ij]>0?-mu*p1+K*s1:-mu*p2+K*s2;
					t2=v[ij]-1>0?mu*u1+tau1:mu*u2+tau2;
					t3=v[ij]+1>0?-mu*u1+tau1:-mu*u2+tau2;
					t4=v[ij]-Kmf>0?-mu*Kmf*v1+p1+s1:-mu*Kmf*v2+p2+s2;
					t5=v[ij]+Kmf>0?mu*Kmf*v1+p1+s1:mu*Kmf*v2+p2+s2;

					uy=0.5*mu_inv*(t2-t3);
					vy=0.5*mu_inv/Kmf*(t5-t4);
					py=(-t1+0.5*K*(t4+t5))/(K+mu);
					sy=(t1+0.5*mu*(t4+t5))/(K+mu);
					tauy=0.5*(t2+t3);
#else
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
#endif
#ifdef FULL_ENO
					cu[ij]=-u[ij]*ux-v[ij]*uy+(-px+sx+tauy)*mu_inv;
					cv[ij]=-u[ij]*vx-v[ij]*vy+(-py-sy+taux)*mu_inv;
					cp[ij]=-u[ij]*px-v[ij]*py-K*(ux+vy);
#else
					cu[ij]=-u[ij]*u1-v[ij]*u2+(-px+sx+tauy)*mu_inv;
					cv[ij]=-u[ij]*v1-v[ij]*v2+(-py-sy+taux)*mu_inv;
					cp[ij]=-u[ij]*p1-v[ij]*p2-K*(ux+vy);
#endif
					magr=sqrt(magn=s[ij]*s[ij]+tau[ij]*tau[ij]);
#ifdef NEW_STZ_MODEL
					cha=dt*nu*exp(-1/chi[ij])*qs_new(magr,chi[ij]);
#else
					cha=dt*nu*exp(-1/chi[ij])*qs(magr);
#endif
#ifndef IMPLICIT
					if(cha>0.45) {
						cout << "Truncate cha " << ij%m << " " << ij/m << endl;
						cha=0.45;
					}
#ifdef FULL_ENO
					cchi[ij]=-u[ij]*chix-v[ij]*chiy+cha*(chi_inf-chi[ij])*c_0_inv;
					magr=mu*cha/(fabs(magr)>1e-8?magn:1);
					cs[ij]=-u[ij]*sx-v[ij]*sy-tomega*tau[ij]-s[ij]*magr+mu*(ux-vy);
					ctau[ij]=-u[ij]*taux-v[ij]*tauy+tomega*s[ij]-tau[ij]*magr+mu*(uy+vx);
#else
					cchi[ij]=-u[ij]*chi1-v[ij]*chi2+cha*(chi_inf-chi[ij])*c_0_inv;
					magr=mu*cha/(fabs(magr)>1e-8?magn:1);
					cs[ij]=-u[ij]*s1-v[ij]*s2-tomega*tau[ij]-s[ij]*magr+mu*(ux-vy);
					ctau[ij]=-u[ij]*tau1-v[ij]*tau2+tomega*s[ij]-tau[ij]*magr+mu*(uy+vx);
#endif
#else
#ifdef FULL_ENO
					cchi[ij]=(chi[ij]-u[ij]*chix-v[ij]*chiy+cha*chi_inf*c_0_inv)/(1+cha*c_0_inv);
					magr=1/(1+mu*cha/(fabs(magr)>1e-8?magn:1));
					cs[ij]=(s[ij]-u[ij]*sx-v[ij]*sy-tomega*tau[ij]+mu*(ux-vy))*magr;
					ctau[ij]=(tau[ij]-u[ij]*taux-v[ij]*tauy+tomega*s[ij]+mu*(uy+vx))*magr;
#else
					cchi[ij]=(chi[ij]-u[ij]*chi1-v[ij]*chi2+cha*chi_inf*c_0_inv)/(1+cha*c_0_inv);
					magr=1/(1+mu*cha/(fabs(magr)>1e-8?magn:1));
					cs[ij]=(s[ij]-u[ij]*s1-v[ij]*s2-tomega*tau[ij]+mu*(ux-vy))*magr;
					ctau[ij]=(tau[ij]-u[ij]*tau1-v[ij]*tau2+tomega*s[ij]+mu*(uy+vx))*magr;
#endif
#endif
#ifdef VISCOSITY
					cu[ij]+=visc*(uxx+hyy*(uu-2*u[ij]+ud));
					cv[ij]+=visc*(vxx+hyy*(vu-2*v[ij]+vd));
#endif
				}
			}

			//Apply the changes
			for(ij=j=0;j<n;j++) for(i=min_i;i<=max_i;i++) {
				ij=i+m*j;
				if(c[ij]<4) {
					u[ij]+=cu[ij];cu[ij]=0;
					v[ij]+=cv[ij];cv[ij]=0;
					p[ij]+=cp[ij];cp[ij]=0;
#ifndef IMPLICIT
					s[ij]+=cs[ij];cs[ij]=0;
					tau[ij]+=ctau[ij];ctau[ij]=0;
					chi[ij]+=cchi[ij];cchi[ij]=0;
#else
					s[ij]=cs[ij];cs[ij]=0;
					tau[ij]=ctau[ij];ctau[ij]=0;
					chi[ij]=cchi[ij];cchi[ij]=0;
#endif
				}
			}

			nucleate_void();
//			cout << wallu << " " << wallx << " " << min_i << " " << max_i << endl;

		//	cout << "Time: " << time << endl;
			print_maxima();
		}

//		cout << time << endl;

		cout << "#Output frame " << k << endl;

		//Clean up grid before output
		write_files(k);
	}
	wallfile.close();
}

/** Computes the magnitudes of the deformation rate tensors.
void sbsim::def_rate_computation() {
	const double small_dev_cutoff=1e-12;
	int i,j,ij;
	double ux,vx,uy,vy,chit,t1,t2,t3,t4;

	// Carry out extrapolation of u and v fields
	ls.extrapolate_fields();

	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) {
		if(c[ij]<4) {

			// Compute centered and one-sided derivatives
			if(i>0) {
				if(i<m-1) {ux=0.5*xsp*(u[ij-1]-u[ij+1]);vx=0.5*xsp*(v[ij-1]-v[ij+1]);}
				else {ux=xsp*(u[ij-1]-u[ij]);vx=xsp*(v[ij-1]-v[ij]);}
			} else {ux=xsp*(u[ij]-u[ij+1]);vx=xsp*(v[ij]-v[ij+1]);}
			if(j>0) {
				if(j<n-1) {uy=0.5*ysp*(u[ij-m]-u[ij+m]);vy=0.5*ysp*(v[ij-m]-v[ij+m]);}
				else {uy=ysp*(u[ij-m]-u[ij]);vy=ysp*(v[ij-m]-v[ij]);}
			} else {uy=ysp*(u[ij]-u[ij+m]);vy=ysp*(v[ij]-v[ij+m]);}

			// Total deformation rate

			cp[ij]=sqrt(0.5*(ux*ux+vy*vy)+0.25*(uy+vx)*(uy+vx));

			// Plastic deformation rate
			magr=sqrt(magn=s[ij]*s[ij]+tau[ij]*tau[ij]);
			cha=nu*exp(-1/chi[ij])*qs(magr);
			magr=1/(1+mu*cha/(fabs(magr)>1e-8?magn:1));
			t4=-log(adaptive_plastic_term(ij,chit,dt))*idt;
			cs[ij]=(s[ij]-u[ij]*s1-v[ij]*s2+tomega*tau[ij]+mu*(ux-vy))*magr;
			ctau[ij]=(tau[ij]-u[ij]*tau1-v[ij]*tau2-tomega*s[ij]+mu*(uy+vx))*magr;
			cs[ij]=t4;

			// Elastic deformation rate
			t3=dev(ij);
			t4/=(t3<small_dev_cutoff?small_dev_cutoff:t3);
			t1=ux-s[ij]*t4;
			t2=0.5*(uy+vx)-tau[ij]*t4;
			t3=vy+s[ij]*t4;
			ctau[ij]=sqrt(0.5*(t1*t1+t3*t3)+t2*t2);
		} else cp[ij]=cs[ij]=ctau[ij]=0;
	}
}*/

void sbsim::nucleate_void() {
#ifdef FULL_ENO
	const double pressure_cutoff=-3.0;
#else
	const double pressure_cutoff=-4.0;
#endif
	const double void_size=0.9*dx;
	bool any_voids=false;
	int i,j,ij;
	for(ij=j=0;j<n;j++) for(i=0;i<m;i++,ij++) if(c[ij]<4&&p[ij]<pressure_cutoff) {
		int i2,j2,ij2;double xs,ys,phit;
		for(ij2=j2=0;j2<n;j2++) for(i2=0;i2<m;i2++,ij2++) {
			xs=dx*(i2-i);ys=dy*(j2-j);
			phit=void_size-sqrt(xs*xs+ys*ys);
			if(phit>phi[ij2]) phi[ij2]=phit;
			any_voids=true;
		}
	}
	if(any_voids) {
		cout << "#Void introduced" << endl;
		ls.build_band();
	}
}

void sbsim::print_maxima() {
	double minu,maxu,minp,maxp,minchi,maxchi;
	int ij=0;
	while(c[ij]>=4) ij++;
	minu=maxu=u[ij];
	minp=maxp=p[ij];
	minchi=maxchi=chi[ij];
	for(;ij<mn;ij++) if(c[ij]<4) {
		if(u[ij]<minu) minu=u[ij];
		if(u[ij]>maxu) maxu=u[ij];
		if(p[ij]<minp) minp=p[ij];
		if(p[ij]>maxp) maxp=p[ij];
		if(chi[ij]<minchi) minchi=chi[ij];
		if(chi[ij]>maxchi) maxchi=chi[ij];
	}
	cout << min_i << " " << max_i << " ";
	printf("Ext: %f %f %f %f %f %f %f\n",time,minu,maxu,minp,maxp,minchi,maxchi);
}

void sbsim::simple_extrapolation(int i,int i2,int i3) {
	int jm;
	for(jm=0;jm<mn;jm+=m) {
	//	u[jm+i]=2*u[jm+i2]-u[jm+i3];
	//	v[jm+i]=2*v[jm+i2]-v[jm+i3];
		p[jm+i]=2*p[jm+i2]-p[jm+i3];
		s[jm+i]=2*s[jm+i2]-s[jm+i3];
		tau[jm+i]=2*tau[jm+i2]-tau[jm+i3];
		chi[jm+i]=2*chi[jm+i2]-chi[jm+i3];
	}
}

inline double sbsim::eno2(double p0,double p1,double p2,double p3) {
	return abs(p0-2*p1+p2)>abs(p1-2*p2+p3)?3*p1-4*p2+p3:p0-p2;
}

inline double sbsim::qs(double s) {
	return fabs(s)>1?(fabs(s)-1)*(fabs(s)-1):0;
}

inline double sbsim::qs_new(double s,double cchi) {
	s=fabs(s);
	double temp=tanh(s/cchi)-1/s;
	return temp>0?cosh(s/cchi)*temp:0;
}

bool sbsim::compute_bc_deriv(double &uu,double &vu,double &pu,double &su,double &tauu,double &chiu,int i,int j,int ni,int nj) {
	if (ni<min_i||ni>max_i) {
		int ij=i+m*j,bij=2*i-ni+m*j;int hm=m>>1;
		double z=wallx-abs(ax+dx*i);
//		if(z>0.75) {
//			uu=((i>hm?wallu:-wallu)-u[ij]*(1-z))/z;
//			vu=((i>hm?wallv:-wallv)-v[ij]*(1-z))/z;
//		} else {
			uu=(2*(i>hm?wallu:-wallu)-u[bij]*(1-z))/(1+z);
//		}
#ifdef SLIDING_BOUNDARY
		vu=(2*v[ij]-v[bij]);
		tauu=-tau[bij]*(1-z)/(1+z);
#else
		vu=(2*(i>hm?wallv:-wallv)-v[bij]*(1-z))/(1+z);
		tauu=2*tau[ij]-tau[bij];
#endif
		pu=2*p[ij]-p[bij];
		su=2*s[ij]-s[bij];
		chiu=2*chi[ij]-chi[bij];
		return false;
	} else if (nj<0||nj>=n) {
		fatal_error("No boundary condition set for top and bottom boundaries",1);
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
		//		if(z>0.75) {
		//			pu=(pp-p[ij]*(1-z))/z;
		//			su=(sp-s[ij]*(1-z))/z;
		//			tauu=(taup-tau[ij]*(1-z))/z;
		//		} else {
					int bi=2*i-ni,bj=2*j-nj;
					int bij=bi+m*bj;
					pu=(2*pp-p[bij]*(1-z))/(1+z);
					su=(2*sp-s[bij]*(1-z))/(1+z);
					tauu=(2*taup-tau[bij]*(1-z))/(1+z);
		//		}
			} else {
				pu=p[nij];su=s[nij];tauu=tau[nij];
			}
			return false;
		}
	}
}

inline double sbsim::phi_x(int i,int ij) {
	return (i==0?phi[ij+1]-phi[ij]:(i==m-1?phi[ij]-phi[ij-1]:(phi[ij+1]-phi[ij-1])*0.5))*xsp;
}

inline double sbsim::phi_y(int i,int ij) {
	return (ij<m?phi[ij+m]-phi[ij]:(ij>=mn-m?phi[ij]-phi[ij-m]:(phi[ij+m]-phi[ij-m])*0.5))*ysp;
}

void sbsim::clean() {
	for(int ij=0;ij<mn;ij++) {
		if(c[ij]>3) {
			u[ij]=v[ij]=p[ij]=s[ij]=tau[ij]=chi[ij]=0;
		} else {
			if (ij%m<min_i) {
				u[ij]=-wallu;
				v[ij]=-wallv;
				p[ij]=s[ij]=tau[ij]=0;chi[ij]=0.074;
			} else if (ij%m>max_i) {
				u[ij]=wallu;
				v[ij]=wallv;
				p[ij]=s[ij]=tau[ij]=0;chi[ij]=0.074;
			}
		}
	}
}

void sbsim::write_files(int k) {
//	clean();
	wallfile << k << " " << wallx << endl;
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
}

template<class T>
void sbsim::output(const char* filename,const char *suffix,T *array,const int sn) {
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

double sbsim::velocity(double x,double y) {
	double x0,y0;
	ls.bicubic(x,y,x0,y0);
	if(abs(x)>wallx) {
		double ut(x>0?wallu:-wallu);
		double temp=-ut*x0-ls.bicubic(v,x>0?wallx:-wallx,y)*y0;
		return temp;
	}/*
#ifdef SLIDING_BOUNDARY
//		x=x>0?wallx-dx*0.5:dx*0.5-wallx;
		x=x>0?wallx:-wallx;
		double x
#else
		return 0;
#endif
//	{
//		return x>0?-vert_speed*y0-horiz_speed*x0:vert_speed*y0+horiz_speed*x0;
//	}*/
	double temp=-ls.bicubic(u,x,y)*x0-ls.bicubic(v,x,y)*y0;
	return temp;
}

void sbsim::fatal_error(const char *p,int code) {
	cerr << p << endl;
	cerr << "Current time: " << time << endl;
	exit(code);
}
