#include <cstdio>
#include <cstdlib>

#include "fileinfo.hh"

class instab_test : public fileinfo {
	public:
		FILE *fp;
		double time;
		double sbar;
		double chiv;
		const char *filename;
		instab_test(const char *filename_) : fileinfo(filename_),
       			filename(filename_) {}
		void solve();
	private:
		void output_line() {
			fprintf(fp,"%.14g %.14g %.14g %.14g %.14g %.14g\n",
				time,time*t_scale,typical_ki*1e-6*(lambda+gamma*time),
				sbar,chiv,chiv*TZ);
		}
		double adaptive_plastic_term(double sbar,double &chiv,double dt);
};

void instab_test::solve() {
	const double small_dev_cutoff=1e-12;
	const double dt=0.001/t_scale;
	const int steps=100;
	const int nevery=100;
	int i,j;
	double dplas,dchi1,dchi2;
	double kfac=(typical_ki*gamma)/sqrt(varrho)/yield_stress;
	char buf[256];

	sprintf(buf,"%s.ist",filename);
	fp=fopen(buf,"w");
	if(fp==NULL) fputs("Can't open output file\n",stderr);

	time=0;sbar=0;chiv=chi0;
	output_line();
	for(i=0;i<steps;i++) {
		for(j=0;j<nevery;j++) {
			dplas=2*mu*t_scale*stz->Dplastic(sbar,chiv,dchi1,dchi2)
				/(sbar<small_dev_cutoff?small_dev_cutoff:sbar);

			sbar+=dt*kfac;
			sbar/=1+dt*dplas;

			chiv=(chiv+(dchi1*chi_inf+dchi2*theta)*t_scale*dt)/(1+dt*t_scale*(dchi1+dchi2));
		}
		time+=nevery*dt;
		output_line();
	}
	fclose(fp);
}

double instab_test::adaptive_plastic_term(double sbar,double &chiv,double dt) {
	const double small_dev_cutoff=1e-12;
	bool adapt=true;
	double dtl=dt,osbar=sbar,dplas,dchi1,dchi2,adt;
	if(sbar<small_dev_cutoff) return 1;
	do {
		dplas=2*t_scale*mu*stz->Dplastic(sbar,chiv,dchi1,dchi2)
		      /(sbar<small_dev_cutoff?small_dev_cutoff:sbar);

		if(abs(dtl*dplas)>adapt_fac) {
			adt=adapt_fac/abs(dplas);
			dtl-=adt;
		} else {
			adapt=false;
			adt=dtl;
		}

		sbar/=1+adt*dplas;
		if(!no_chidot) chiv=(chiv+(dchi1*chi_inf+dchi2*theta)*t_scale*adt)/(1+adt*t_scale*(dchi1+dchi2));
	} while(adapt);

	return sbar/osbar;
}

int main(int argc,char **argv) {
	if(argc!=2) {
		fprintf(stderr,"Usage: ./instab_test <sim_name>\n");
		return 1;
	}

	instab_test ist(argv[1]);
	ist.solve();
}
