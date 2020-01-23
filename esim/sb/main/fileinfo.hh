#ifndef FILEINFO_HH
#define FILEINFO_HH

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "stz_model.hh"

class fileinfo {
	public:
		static const int buf_size=512;
		char buf[buf_size];
		const char *filename;
		double chi0;
		double random_chi0_pct;
		double random_chi0_len;
		double T;
		double TZ;
		double sim_size;
		double output_size;
		int m;
		unsigned int fflags;
		double adapt_fac;
		double displacement;
		double viscosity;
		double filter_stress;
		double varrho;
		double youngs_modulus;
		double poisson;
		double lambda;
		double gamma;
		double typical_ki;
		double rho0;
		double yield_stress;
		double init_void_pos;
		double init_void_rad;
		double void_nucl_p;
		double void_nucl_rad;
		double void_nucl_tfac;
		double void_nucl_rand;
		double sec_size;
		double c0;
		double tau0;
		double kappa;
		double Delta;
		double Omega;
		double eps0;
		double chi_inf;
		double chi_len;
		double tmult;
		double max_ki;
		int frames;
		double frho_a;
		double alpha0;
		double T0;
		double TA;
		double T1;
		double theta;
		double trunc_modify_param;
		double bc_smooth;
		double mixity;
		double mix_i,mix_ii;
		int a_substep;
		int ntrace;
		int truncate;
		int trunc_modify;
		bool athermal;
		bool nonlinear;
		bool init_void;
		bool void_nucl;
		bool random_chi0;
		bool output_section;
		bool no_chidot;
		bool limit_output;
		bool irwin_no_fac2;
		bool adaptive_ts;
		double K;
		double mu;
		double Kappa;
		double Kinit;
		double Factor;
		double Omegaeps0;
		double rho;
		double t_scale;
		stz_dynamics *stz;
		fileinfo(const char *filename_);
		~fileinfo();
		void digest();
		void output_html(FILE *fp=stdout);
	protected:
		/** Tests to see if two strings are equal.
		 * \param[in] p1 a pointer to the first string.
		 * \param[in] p2 a pointer to the second string.
		 * \return True if they are equal, false otherwise. */
		inline bool se(const char *p1,const char *p2) {
			return strcmp(p1,p2)==0;
		}
	private:
		void calculate_constants();
		/** Finds the next token in a string and interprets it as a
		 * double precision floating point number. If none is availble,
		 * it gives an error message.
		 * \param[in] ln the current line number. */
		inline double next_double(int ln) {
			return atof(next_token(ln));
		}
		/** Finds the next token in a string, interprets it as a double
		 * precision floating point number, and checks that there are
		 * no subsequent values.
		 * \param[in] ln the current line number. */
		inline double final_double(int ln) {
			double temp=next_double(ln);
			check_no_more(ln);
			return temp;
		}
		char* next_token(int ln);
		void check_no_more(int ln);
		void check_invalid(double val,const char *p);
};

#endif
