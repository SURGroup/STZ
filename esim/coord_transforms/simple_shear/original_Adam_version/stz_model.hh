#ifndef STZ_MODEL_HH
#define STZ_MODEL_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

class stz_dynamics {
    public:
        const double TZ;
        const double chi_inf;
        stz_dynamics(double TZ_,double chi_inf_) : TZ(TZ_), chi_inf(chi_inf_) {}
        virtual ~stz_dynamics() {}
        virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) = 0;
};

class stz_dynamics_chr : public stz_dynamics {
    public:
        const double c0;
        const double tau0;
        const double kappa;
        const double Delta;
        const double Omegaeps0;
        const double thetainv;
        const double rho;
        const double s0;
        stz_dynamics_chr(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
                double Omegaeps0_,double chi_inf_,double theta,double rho_) :
            stz_dynamics(TZ_,chi_inf_), c0(c0_), tau0(tau0_),
            kappa(kappa_), Delta(Delta_), Omegaeps0(Omegaeps0_),
            thetainv(1/theta), rho(rho_), s0(1) {}
        virtual ~stz_dynamics_chr() {}
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

class stz_dynamics_adam_original : public stz_dynamics_chr {
        public:
        const double mu;
        const double mu_inv;
                stz_dynamics_adam_original(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
                        double Omegaeps0_,double chi_inf_,double theta,double rho_,double mu_) :
                        stz_dynamics_chr(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_),
                   mu(mu_), mu_inv(1/mu) {};
                ~stz_dynamics_adam_original() {};
                virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) {

          //compute Dplastic
          double bfac=exp(-1/chi),
             Dpl=bfac*exp(1/chi_inf)*sbar*mu_inv;

                        // Calculate my plasticity model
                        dchi2=0;
            dchi1=2*tau0*Dpl*sbar/(s0*c0);
                        return Dpl;
                }
};

class stz_dynamics_linear : public stz_dynamics_chr {
    public:
        int t_model;
        stz_dynamics_linear(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
            double Omegaeps0_,double chi_inf_,double theta,double rho_,int t_model_=-1) :
            stz_dynamics_chr(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_),
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

class stz_dynamics_nonlinear : public stz_dynamics_chr {
    public:
        stz_dynamics_nonlinear(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
            double Omegaeps0_,double chi_inf_,double theta,double rho_) :
            stz_dynamics_chr(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_) {};
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

class stz_dynamics_linear_athermal : public stz_dynamics_chr {
    public:
        int t_model;
        stz_dynamics_linear_athermal(double TZ_,double c0_,double tau0_,double kappa_,double Delta_,
            double Omegaeps0_,double chi_inf_,double theta,double rho_,int t_model_=-1) :
            stz_dynamics_chr(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_),
            t_model(t_model_) {};
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
                    sC=0.5*Omegaeps0*sbar/Delta;
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
};

class stz_dynamics_nonlinear_athermal : public stz_dynamics_chr {
    public:
        stz_dynamics_nonlinear_athermal(double TZ_,double c0_,double tau0_,double kappa_,
                double Delta_,double Omegaeps0_,double chi_inf_,
                double theta,double rho_) :
            stz_dynamics_chr(TZ_,c0_,tau0_,kappa_,Delta_,Omegaeps0_,chi_inf_,theta,rho_) {};
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

class stz_dynamics_adam : public stz_dynamics {
    public:
        const double tau0;
        const double ep;
        const double c0;
        stz_dynamics_adam(double TZ_,double chi_inf_,double tau0_,double ep_,double c0_) :
                    stz_dynamics(TZ_,chi_inf_), tau0(tau0_), ep(ep_), c0(c0_) {};
            ~stz_dynamics_adam() {};
            virtual double Dplastic(double sbar,double chi,double &dchi1,double &dchi2) {

        // If the material is below yield, then there is zero
        // deformation
        if(sbar<1) {
            dchi1=dchi2=0;
            return 0;
        }

        // Compute Dplastic - no need to do abs(sbar) since
        // sbar is always positive
        double C=-2+sbar+exp(-sbar)*(2+sbar),
               Dpl=(ep/tau0)*C*(1-1/sbar)*exp(-1/chi);

        // Compute changes to chi
        dchi1=2*Dpl*sbar/c0;
        dchi2=0;
        return Dpl;
    }
};

#endif
