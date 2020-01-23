#ifndef STZ_MODEL_HH
#define STZ_MODEL_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

class stz_dynamics {
	public:
		const double rt3;
		const double TZ;
		const double c0;
		const double tau0;
		const double kappa;
		const double Delta;
		const double Omegaeps0;
		const double chi_inf;
		const double thetainv;
		const double rho;
		const double s0;
		stz_dynamics(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
				double Omegaeps0_,double chi_inf_,double theta,double rho_) :
			rt3(sqrt(3.0)), TZ(TZ_), c0(c0_), tau0(tau0_),
			kappa(kappa_), Delta(Delta_), Omegaeps0(Omegaeps0_),
			chi_inf(chi_inf_), thetainv(1/theta), rho(rho_), s0(1) {}
		virtual ~stz_dynamics() {}
		virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) = 0;
	protected:
		inline double Dplastic_common(double sbar,double chi,double &dchi1,double &dchi2,
				double sC1,double sC2,double sC,double sT) {

			// Compute M function, Dplastic and Gamma
			double sbart=sbar>1e-12?sbar:1e-12,bfac=exp(-1/chi),
			       pM=1+sbart*sT/s0+rho/(2*sC2),
			       M=s0/(2*sbart)*(pM-sqrt(pM*pM-4*sbart*sT/s0)),
			       Dpl=bfac*sC*(sT-M)/tau0,
			       Gamma=2*tau0*Dpl*sbar/(s0*bfac);

			// Compute changes to chi
			bfac/=tau0*c0;
			dchi1=bfac*Gamma;
			dchi2=bfac*kappa*sC1*rho;

			// Return the plastic deformation
			return Dpl;
		}
};

class stz_dynamics_linear : public stz_dynamics {
	public:
		int t_model;
		stz_dynamics_linear(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
			double Omegaeps0_,double chi_inf_,double theta,double rho_,int t_model_=-1) :
			stz_dynamics(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_),
       			t_model(t_model_) {};
		~stz_dynamics_linear() {};
		virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) {

			// Calculate STZ-related functions
			double sC1=exp(-Delta*thetainv);
			double sC2=cosh(Omegaeps0*sbar*thetainv),sC=sC1*sC2,
			       sT=tanh(Omegaeps0*sbar*thetainv);

			// Call common computation routine
			return Dplastic_common(sbar,chi,dchi1,dchi2,sC1,sC2,sC,sT);
		}
};

class stz_dynamics_nonlinear : public stz_dynamics {
	public:
		stz_dynamics_nonlinear(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
			double Omegaeps0_,double chi_inf_,double theta,double rho_) :
			stz_dynamics(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_) {};
		~stz_dynamics_nonlinear() {};
		virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) {

			// Calculate STZ-related functions
			double x=Omegaeps0*sbar/Delta,
			       sC1=exp(-Delta*thetainv),
			       sC=exp(-Delta*cosh(x)*thetainv)*cosh(Delta*sinh(x)*thetainv),
			       sC2=sC/sC1,
			       sT=tanh(Delta*sinh(x)*thetainv);

			// Call common computation routine
			return Dplastic_common(sbar,chi,dchi1,dchi2,sC1,sC2,sC,sT);
		}
};

class stz_dynamics_linear_athermal : public stz_dynamics {
	public:
		/** The trunction model to apply. */
		int t_model;
		/** For truncation model 3, the intercept of the linear truncation function. */
		double t_icpt;
		/** For truncation model 3, the gradient of the linear truncation function. */
		double t_grad;
		stz_dynamics_linear_athermal(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
			double Omegaeps0_,double chi_inf_,double theta,double rho_,int t_model_=-1) :
			stz_dynamics(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_),
			t_model(t_model_), t_icpt(0), t_grad(0.5*Omegaeps0/Delta) {};
		~stz_dynamics_linear_athermal() {};
		virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) {
			double tmp,sC,sT;

			// There is no increase in chi due to thermal effects
			dchi2=0;

			// If the material is below yield, then there is zero
			// deformation
			if(sbar<1) {
				dchi1=0;
				return 0;
			}

			// Calculate STZ-related functions, applying the
			// truncation if required
			if(t_model==-1||Delta>=Omegaeps0*sbar) sC=0.5*(exp((Omegaeps0*sbar-Delta)*thetainv)+exp((-Omegaeps0*sbar-Delta)*thetainv));
			else switch(t_model) {
				case 0:
					sC=0.5*(exp((-Omegaeps0*sbar-Delta)*thetainv)+1);
					break;
				case 1:
					tmp=Omegaeps0*sbar/Delta;
					sC=0.5*(exp((-Omegaeps0*sbar-Delta)*thetainv)+sqrt(0.5*(1+tmp*tmp)));
					break;
				case 2:
					tmp=Omegaeps0*sbar/Delta;
					sC=0.5*(exp((-Omegaeps0*sbar-Delta)*thetainv)+(0.5*(1+tmp*tmp)));
					break;
				case 3:
					sC=t_icpt+sbar*t_grad;
					break;
				case 4:
					sC=0.5*(exp((-Omegaeps0*sbar-Delta)*thetainv)+(Omegaeps0*sbar-Delta)*thetainv+1);break;
				default:
					sC=0;
			}
			sT=t_model!=3?tanh(Omegaeps0*sbar*thetainv):1;

			// Compute Dplastic and Gamma
			double bfac=exp(-1/chi),
			       Dpl=bfac*sC*(sT-s0/sbar)/tau0,
			       Gamma=2*tau0*Dpl*sbar/(s0*bfac);

			// Compute changes to chi
			dchi1=bfac*Gamma/(tau0*c0);

			// Return the plastic deformation
			return Dpl;
		}
		/** Modifies the slope of truncation model 3, scaling it by a
		 * multiple from the standard value.
		 * \param[in] p the multiplicative factor. */
		void trunc_modify_mult(double p) {
			t_grad=0.5*Omegaeps0/Delta*p;set_icpt();
		}
		/** Modifies the slope of truncation model 3, scaling it
		 * between the standard slope (for p=0) and the slope for which
		 * it would be differentiable at the crossover (p=1).
		 * \param[in] p the scaling factor. */
		void trunc_modify_scale(double p) {
			t_grad=(1-p)*(0.5*Omegaeps0/Delta)+p*(0.5*(1-exp(-2*Delta*thetainv))*Omegaeps0*thetainv);
			set_icpt();
		}
	private:
		inline void set_icpt() {
			t_icpt=0.5*(1+exp(-2*Delta*thetainv))-Delta/Omegaeps0*t_grad;
		}
};

class stz_dynamics_nonlinear_athermal : public stz_dynamics {
	public:
		stz_dynamics_nonlinear_athermal(double TZ_,double c0_,double tau0_,double kappa_,
				double Delta_,double Omegaeps0_,double chi_inf_,
				double theta,double rho_) :
			stz_dynamics(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_) {};
		~stz_dynamics_nonlinear_athermal() {};
		virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) {

			// There is no increase in chi due to thermal effects
			dchi2=0;

			// If the material is below yield, then there is zero
			// deformation
			if(sbar<1) {
				dchi1=0;
				return 0;
			}

			// Calculate STZ-related functions
			double x=Omegaeps0*sbar/Delta,
			       sC=0.5*(exp(-Delta*thetainv*exp(x))+exp(-Delta*thetainv*exp(-x))),
			       sT=tanh(Delta*sinh(x)*thetainv);

			// Compute Dplastic and Gamma
			double bfac=exp(-1/chi),
			       Dpl=bfac*sC*(sT-s0/sbar)/tau0,
			       Gamma=2*tau0*Dpl*sbar/(s0*bfac);

			// Compute changes to chi
			dchi1=bfac*Gamma/(tau0*c0);

			// Return the plastic deformation
			return Dpl;
		}
};

class stz_dynamics_damage {
	public:
		const double alpha1;
		const double alpha2;
		const double rho0;
		const double k0;
		const double gamma_g;
		const double pfac;
		const double chifac;
		stz_dynamics_damage(double alpha1_,double alpha2_,double rho0_,double k0_,
				double E,double gamma_g_,double V0,double C,double aa_init) :
				alpha1(alpha1_), alpha2(alpha2_), rho0(rho0_),
				k0(k0_), gamma_g(gamma_g_),
				pfac(2*E*gamma_g/3.1415926535897932384626433832795),
				chifac(C/(V0*aa_init*aa_init*aa_init)) {}
		double Dplastic(double sbar,double chi,double p,double aa,double &dchi1,double &daa) {
			double itau,s0,Dpl,pth,epthp,c0;

			// Yield stress
			s0=alpha1+alpha2*p;

			// Skip if no yielding
			if(sbar<s0) {
				dchi1=daa=0;
				return 0;
			}

			// Inverse timescale
			itau=sqrt(p/rho0)/aa;

			// Plastic deformation
			Dpl=itau*exp((-1+sbar/p)/chi)*(1-s0/sbar);

			// P_th calculations
			pth=sqrt(pfac/aa);
			epthp=exp(-pth/p);

			// Change in grain size
			daa=-k0*epthp*sbar*Dpl/gamma_g*aa*aa;

			// Change in effective temperature
			c0=chifac*aa*aa*aa;
			dchi1=sbar*Dpl/p*(1-k0*epthp)/c0;
			if(rand()%10000==0) printf("YO %g %g %g %g %g Dpl=%g --- %g %g --- %g %g\n",sbar,p,s0,exp((-1+sbar/p)/chi),dchi1,Dpl,aa,daa,pth,1-k0*epthp);

			return Dpl;
		}
};

#endif
