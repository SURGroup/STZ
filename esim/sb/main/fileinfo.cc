#include "fileinfo.hh"

/** The fileinfo constructor parses a simulation setup file, and sets all of the
 * required constants. It calculates dependent constants and creates the STZ model
 * class.
 * \param[in] fn the filename of the setup file. */
fileinfo::fileinfo(const char *filename_) : filename(filename_), chi0(-1),
	T(-1), TZ(-1), sim_size(-1), m(-1), fflags(0), adapt_fac(-1),
	displacement(-1), viscosity(-1), filter_stress(0), varrho(-1),
	youngs_modulus(-1), poisson(-1), lambda(-1), gamma(-1), typical_ki(-1),
	rho0(-1), yield_stress(-1), init_void_pos(-1), init_void_rad(-1),
	void_nucl_p(-1), void_nucl_rad(-1), void_nucl_tfac(-1), c0(-1),
	tau0(-1), kappa(-1), Delta(-1), Omega(-1), eps0(-1), chi_inf(-1),
	chi_len(0), tmult(-1), max_ki(-1), frames(0), frho_a(-1), alpha0(-1),
	T0(-1), TA(-1), T1(-1), theta(-1), bc_smooth(0), mixity(0),
	mix_i(1), mix_ii(0), a_substep(1), truncate(-1), trunc_modify(-1),
	athermal(false), nonlinear(false), init_void(false),
	void_nucl(false), random_chi0(false), output_section(false),
	no_chidot(false), limit_output(false), irwin_no_fac2(false),
	adaptive_ts(false) {
	int ln=1;

	// Append ".fin" to filename and open file
	sprintf(buf,"%s.fin",filename);
	FILE *f=fopen(buf,"r");
	if(f==NULL) {
		fprintf(stderr,"Can't open input file '%s'\n",buf);
		exit(1);
	}
	char *bp;

	// Read in lines of the input file
	while(!feof(f)) {
		if(fgets(buf,buf_size,f)==NULL) break;

		// Locate comments and remove by replacing comment character
		// with a null character
		bp=buf;
		while((*bp)!=0) {
			if(*bp=='#') {*bp=0;break;}
			bp++;
		}

		// Search for a keyword, and skip if none found
		bp=strtok(buf," \t\n");
		if(bp==NULL) continue;

		if(se(bp,"chi0")) chi0=final_double(ln);
		else if(se(bp,"random_chi0")) {
			random_chi0=true;
			random_chi0_pct=next_double(ln);
			random_chi0_len=final_double(ln);
		}
		else if(se(bp,"T")) T=final_double(ln);
		else if(se(bp,"TZ")) TZ=final_double(ln);
		else if(se(bp,"model")) {
			bp=next_token(ln);
			while(bp!=NULL) {
				if(se(bp,"athermal")) athermal=true;
				else if (se(bp,"thermal")) athermal=false;
				else if (se(bp,"nonlinear")) nonlinear=true;
				else if (se(bp,"linear")) nonlinear=false;
				else {
					fprintf(stderr,"Can't understand model parameter '%s' at line %d\n",bp,ln);
					exit(1);
				}
				bp=strtok(NULL," \t\n");
			}
		} else if(se(bp,"truncate_model")) {
			bp=next_token(ln);
			truncate=atoi(bp);
			if(truncate<0||truncate>4) {
				fprintf(stderr,"Truncation model must be between 0 and 4 at line %d\n",ln);
				exit(1);
			}
		} else if(se(bp,"truncate_modify")) {
			bp=next_token(ln);
			if(se(bp,"mult")) trunc_modify=0;
			else if(se(bp,"scale")) trunc_modify=1;
			else {
				fputs("Invalid truncation modification model\n",stderr);
				exit(1);
			}
			trunc_modify_param=final_double(ln);
		} else if(se(bp,"sim_size")) sim_size=final_double(ln);
		else if(se(bp,"output_size")) {
			limit_output=true;
			output_size=final_double(ln);
		} else if(se(bp,"gridpoints")) {
			bp=next_token(ln);
			m=atoi(bp)-1;
			if(m<=0) {
				fprintf(stderr,"Gridpoints must be positive at line %d\n",ln);
				exit(1);
			}
		} else if(se(bp,"output_section")) {
			output_section=true;
			sec_size=final_double(ln);
		} else if(se(bp,"output")) {
			fflags=0;

			bp=next_token(ln);
			while(bp!=NULL) {
				if(se(bp,"u")) fflags|=1;
				else if(se(bp,"v")) fflags|=2;
				else if(se(bp,"p")) fflags|=4;
				else if(se(bp,"q")) fflags|=8;
				else if(se(bp,"s")) fflags|=16;
				else if(se(bp,"tau")) fflags|=32;
				else if(se(bp,"chi")) fflags|=64;
				else if(se(bp,"tem")) fflags|=128;
				else if(se(bp,"dev")) fflags|=256;
				else if(se(bp,"X")) fflags|=512;
				else if(se(bp,"Y")) fflags|=1024;
				else if(se(bp,"qs")) fflags|=2048;
				else if(se(bp,"VS")) fflags|=4096;
				else if(se(bp,"DS")) fflags|=8192;
				else if(se(bp,"Dtot")) fflags|=16384;
				else if(se(bp,"Dpl")) fflags|=32768;
				else if(se(bp,"Del")) fflags|=65536;
				else if(se(bp,"cc")) fflags|=1<<17;
				else if(se(bp,"d")) fflags|=1<<18;
				else if(se(bp,"c")) fflags|=1<<19;
				else if(se(bp,"phi")) fflags|=1<<20;
				bp=strtok(NULL," \t\n");
			}
		} else if(se(bp,"adapt_fac")) adapt_fac=final_double(ln);
		else if(se(bp,"displacement")) displacement=final_double(ln);
		else if(se(bp,"viscosity")) viscosity=final_double(ln);
		else if(se(bp,"filter_stress")) filter_stress=final_double(ln);
		else if(se(bp,"notch_radius")) varrho=final_double(ln);
		else if(se(bp,"youngs_modulus")) youngs_modulus=final_double(ln);
		else if(se(bp,"poisson_ratio")) poisson=final_double(ln);
		else if(se(bp,"lambda")) lambda=final_double(ln);
		else if(se(bp,"gamma")) gamma=final_double(ln);
		else if(se(bp,"typical_ki")) typical_ki=final_double(ln);
		else if(se(bp,"rho0")) rho0=final_double(ln);
		else if(se(bp,"yield_stress")) yield_stress=final_double(ln);
		else if(se(bp,"init_void")) {
			init_void=true;
			init_void_pos=next_double(ln);
			init_void_rad=final_double(ln);
		} else if(se(bp,"void_nucl")) {
			void_nucl=true;
			void_nucl_p=next_double(ln);
			void_nucl_rad=next_double(ln);
			void_nucl_tfac=next_double(ln);
			void_nucl_rand=final_double(ln);
		} else if(se(bp,"c0")) c0=final_double(ln);
		else if(se(bp,"tau0")) tau0=final_double(ln);
		else if(se(bp,"kappa")) kappa=final_double(ln);
		else if(se(bp,"Delta")) Delta=final_double(ln);
		else if(se(bp,"Omega")) Omega=final_double(ln);
		else if(se(bp,"eps0")) eps0=final_double(ln);
		else if(se(bp,"chi_inf")) chi_inf=final_double(ln);
		else if(se(bp,"chi_len")) chi_len=final_double(ln);
		else if(se(bp,"bc_smooth")) bc_smooth=final_double(ln);
		else if(se(bp,"tmult")) tmult=final_double(ln);
		else if(se(bp,"max_ki")) max_ki=final_double(ln);
		else if(se(bp,"frames")) {
			bp=next_token(ln);
			frames=atoi(bp);
			if(frames<0) {
				fprintf(stderr,"Number of frames must be non-negative at line %d\n",ln);
				exit(1);
			}
			check_no_more(ln);
		} else if(se(bp,"frho")) {
			frho_a=next_double(ln);
			alpha0=next_double(ln);
			T0=next_double(ln);
			TA=next_double(ln);
			T1=final_double(ln);
		} else if(se(bp,"mode_ii")) mixity=90;
		else if(se(bp,"mixity")) mixity=final_double(ln);
		else if(se(bp,"irwin_no_fac2")) irwin_no_fac2=true;
		else if(se(bp,"adaptive_ts")) adaptive_ts=true;
		else if(se(bp,"no_chidot")) no_chidot=true;
		else if(se(bp,"advect_substep")) {
			bp=next_token(ln);
			a_substep=atoi(bp);
			if(a_substep<1) {
				fprintf(stderr,"Advection substep must be greater than 0 at line %d\n",ln);
				exit(1);
			}
			check_no_more(ln);
		}
		else {
			fprintf(stderr,"Unrecognized keyword '%s' at line %d\n",bp,ln);
			exit(1);
		}
		ln++;
	}
	fclose(f);

	// Check all necessary parameters have been set
	if(m==-1) {
		fputs("Gridpoints not set\n",stderr);
		exit(1);
	}
	if(nonlinear&&truncate!=-1) fputs("Warning: truncate keyword ignored with nonlinear model\n",stderr);
	if(fflags==0) fputs("Warning: no output fields specified\n",stderr);
	check_invalid(chi0,"chi0");
	check_invalid(T,"T");
	check_invalid(TZ,"TZ");
	check_invalid(sim_size,"sim_size");
	check_invalid(adapt_fac,"adapt_fac");
	check_invalid(displacement,"displacement");
	check_invalid(viscosity,"viscosity");
	check_invalid(varrho,"notch_radius");
	check_invalid(youngs_modulus,"youngs_modulus");
	check_invalid(poisson,"poisson_ratio");
	check_invalid(lambda,"lambda");
	check_invalid(gamma,"gamma");
	check_invalid(typical_ki,"typical_ki");
	check_invalid(rho0,"rho0");
	check_invalid(yield_stress,"yield_stress");
	check_invalid(c0,"c0");
	check_invalid(tau0,"tau0");
	check_invalid(kappa,"kappa");
	check_invalid(Delta,"Delta");
	check_invalid(Omega,"Omega");
	check_invalid(eps0,"eps0");
	check_invalid(chi_inf,"chi_inf");
	check_invalid(tmult,"tmult");
	check_invalid(max_ki,"max_ki");
	if(frho_a<=0||alpha0<=0||T0<=0||TA<=0||T1<=0) {
		fputs("rho(T) parameters incorrect or not set\n",stderr);
		exit(1);
	}
	if(trunc_modify!=-1&&(!athermal||nonlinear||truncate!=3)) {
		fputs("Truncation modification not valid for this STZ model\n",stderr);
		exit(1);
	}

	calculate_constants();
}

void fileinfo::calculate_constants() {
	const double pi=3.1415926535897932384626433832795,
		     kB=1.3806503e-23;

	// Rescale temperatures
	chi0/=TZ;chi_inf/=TZ;theta=T/TZ;Delta/=TZ;

	// Elastic modulus computation
	K=youngs_modulus/(3*(1-2*poisson))/yield_stress;
	mu=youngs_modulus/(2*(1+poisson))/yield_stress;

	// Irwin kappa and boundary condition constants
	Kappa=3-4*poisson;
	Kinit=typical_ki*lambda/yield_stress*sqrt(1/(2*pi*varrho));
	Factor=(irwin_no_fac2?1:0.5)*typical_ki*gamma/youngs_modulus*(1+poisson)*sqrt(1/(2*pi*varrho));

	// Scaled Omega*epsilon_0 quantity used in STZ model
	Omegaeps0=Omega*1e-30*eps0*yield_stress/(TZ*kB);

	// f(rho) computation
	rho=T>T0?exp(-T1/(T-T0)*exp(-frho_a*(T-T0)/(TA-T0))-alpha0):0;

	// Timescale computation
	t_scale=varrho*sqrt(rho0/(mu*yield_stress));

	// Create class to handle the STZ dynamics
	stz=athermal?(nonlinear?(stz_dynamics*)
		new stz_dynamics_nonlinear_athermal(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,rho):
			        (stz_dynamics*)
		new stz_dynamics_linear_athermal(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,rho,truncate))
		    :(nonlinear?(stz_dynamics*)
		new stz_dynamics_nonlinear(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,rho):
				(stz_dynamics*)
		new stz_dynamics_linear(TZ,c0,tau0,kappa,Delta,Omegaeps0,chi_inf,theta,rho));

	// Mode mixity computation
	mix_i=cos((pi/180.)*mixity);
	mix_ii=sin((pi/180.)*mixity);

	// Apply truncation modification if enabled
	if(trunc_modify==0)
		((stz_dynamics_linear_athermal*) stz)->trunc_modify_mult(trunc_modify_param);
	else if(trunc_modify==1)
		((stz_dynamics_linear_athermal*) stz)->trunc_modify_scale(trunc_modify_param);
}

/** The class destructor frees the dynamically allocated STZ model. */
fileinfo::~fileinfo() {
	delete stz;
}

/** Prints a digest of all of the run setup. */
void fileinfo::digest() {

	// Print information about model and temperatures
	printf("Simulation name   : %s\n"
	       "STZ Model         : %slinear %sthermal\n"
	       ,filename,nonlinear?"non":"",athermal?"a":"");
	if(!nonlinear) {
		if(truncate==-1) puts("Truncation        : none");
		else {
			printf("Truncation model  : %d\n",truncate);
			if(trunc_modify!=-1)
				printf("Truncation modify : %s\nTrunc. parameter  : %g\n",
				       trunc_modify==0?"mult":"scale",trunc_modify_param);
		}
	}
	printf("Mode mixity       : %g degrees\n\n"
	       "T_Z               : %g K\n"
	       "T                 : %g K (%.5g)\n"
	       "chi0              : %g K (%.5g)\n",mixity,TZ,T,theta,chi0*TZ,chi0);
	if(random_chi0) printf("Random chi0       : %g%%, %g µm (%g)\n",random_chi0_pct,random_chi0_len*varrho*1e6,random_chi0_len);

	// Print computational details
	printf("\nSize              : %g\n"
	       "Gridpoints        : %d by %d\n"
	       "Output fields     :",sim_size,m,m);
	if(fflags&1) printf(" u");
	if(fflags&2) printf(" v");
	if(fflags&4) printf(" p");
	if(fflags&8) printf(" q");
	if(fflags&16) printf(" s");
	if(fflags&32) printf(" tau");
	if(fflags&64) printf(" chi");
	if(fflags&128) printf(" tem");
	if(fflags&256) printf(" dev");
	if(fflags&512) printf(" X");
	if(fflags&1024) printf(" Y");
	if(fflags&2048) printf(" qs");
	if(fflags&4096) printf(" VS");
	if(fflags&8192) printf(" DS");
	if(fflags&16384) printf(" Dtot");
	if(fflags&32768) printf(" Dpl");
	if(fflags&65536) printf(" Del");
	if(fflags&131072) printf(" cc");
	if(fflags&262144) printf(" d");
	if(fflags&524288) printf(" c");
	if(fflags&1048576) printf(" phi");
	printf("\nReference map     : %sused\n"
	       "Adaptivity        : %g\n"
	       "Displacement      : %g\n"
	       "Viscosity         : %g\n\n",fflags&15360?"":"not ",
	       adapt_fac,displacement,viscosity);

	// Print initial void information if present
	if(init_void) printf("Init void pos.    : %g\n"
			     "Init void radius  : %g\n\n",init_void_pos,init_void_rad);

	// Print void nucleation information if present
	if(void_nucl) printf("Void nucl. p      : %.4g GPa (%g)\n"
			     "Void nucl. rad    : %g µm (%g dx)\n"
			     "Void nucl. t. fac.: %g\n\n"
			     "Void nucl. rand   : %g\n",
			     void_nucl_p*yield_stress*1e-9,void_nucl_p,
			     varrho/m*2*sim_size*1e6*void_nucl_rad,
			     void_nucl_rad,void_nucl_tfac,void_nucl_rand);

	// Print mechanical constants
	printf("Notch radius      : %g µm\n"
	       "Young's modulus   : %g GPa\n"
	       "Poisson ratio     : %g\n"
	       "lambda & gamma    : %g & %g\n"
	       "Typical K_I       : %.5g MPa sqrt(m)\n"
	       "Density           : %.5g kg m^{-3}\n"
	       "Yield stress      : %.5g GPa\n"
	       "Bulk modulus      : %.5g GPa (%.5g)\n"
	       "Shear modulus     : %.5g GPa (%.5g)\n"
	       "Time scale        : %.5g s\n\n",
	       varrho*1e6,youngs_modulus*1e-9,poisson,lambda,gamma,typical_ki*1e-6,rho0,
	       yield_stress*1e-9,K*yield_stress*1e-9,K,mu*yield_stress*1e-9,mu,t_scale);

	// Print STZ parameters
	printf("c0                : %g\n"
	       "tau0              : %g s\n"
	       "kappa             : %g\n"
	       "Delta             : %g K (%.5g)\n"
	       "Omega             : %g Å^3\n"
	       "eps0              : %g\n"
	       "chi_inf           : %.5g K (%.5g)\n",c0,tau0,kappa,Delta*TZ,Delta,Omega,eps0,chi_inf*TZ,chi_inf);

	// Print rho for thermal models only
	if(!athermal) printf("\nrho(T)=exp(-T1/(T_-T0)*exp(-a*(T_-T0)/(TA-T0))-alpha0)=%.5g\n"
			     "a=%g alpha0=%g T0=%g TA=%g T1=%g\n",rho,frho_a,alpha0,T0,TA,T1);

	// Print run information
	printf("\nTimestep mult.    : %g\n"
	       "Maximum K_I       : %g MPa sqrt(m) (%.5g)\n"
	       "Frames            : %d\n",tmult,max_ki*1e-6,max_ki/typical_ki,frames);
}

/** Outputs parameters about the setup to take up one line of an HTML table.
 * \param[in] fp the file handle to write to. */
void fileinfo::output_html(FILE *fp) {
	fprintf(fp,"<td><acronym title=\"%sinear\">%sL</acronym>",nonlinear?"Nonl":"L",nonlinear?"N":"");
	if(truncate!=-1) fprintf(fp,"<acronym title=\"Truncation model %d\">T%d</acronym>",truncate,truncate);
	fprintf(fp,"&nbsp;<acronym title=\"%shermal\">%sT</acronym></td>"
		   "<td>%g&nbsp;K</td><td>%g&nbsp;K",
		   athermal?"At":"T",athermal?"A":"",
		   T,chi0*TZ);
	if(random_chi0) fprintf(fp,"&nbsp;(%.3g%%)&nbsp;[%.3g]",random_chi0_pct,random_chi0_len);
	fprintf(fp,"</td><td>%g&nbsp;K</td>"
		   "<td>%d</td><td>%g</td><td>%g</td>"
		   "<td>%g</td><td>%.2g</td><td>%g</td><td>%.4g</td>"
		   "<td>%g</td><td>%g</td><td>%g</td><td>%g&nbsp;K</td><td>%g</td><td>%g</td><td>%g&nbsp;K</td><td>%.3g</td>"
		   "<td>%.2g</td>",
		   TZ,
		   m,sim_size,varrho*1e6,
		   youngs_modulus*1e-9,poisson,yield_stress*1e-9,typical_ki*1e-6*gamma/t_scale,
		   c0,tau0,kappa,Delta*TZ,Omega,eps0,chi_inf*TZ,Delta/Omegaeps0,
		   tmult);
}

/** Checks that there are no subsequent values.
 * \param[in] ln the current line number. */
void fileinfo::check_no_more(int ln) {
	if(strtok(NULL," \t\n")!=NULL) {
		fprintf(stderr,"Too many arguments at input line %d\n",ln);
		exit(1);
	}
}

/** Finds the next token in a string, and if none is availble, gives an error
 * message.
 * \param[in] ln the current line number. */
char* fileinfo::next_token(int ln) {
	char *temp=strtok(NULL," \t\n");
	if(temp==NULL) {
		fprintf(stderr,"Not enough arguments at input line %d\n",ln);
		exit(1);
	}
	return temp;
}

/** Checks that a parameter is a valid positive value.
 * \param[in] val the value to check.
 * \param[in] p the name of the value. */
void fileinfo::check_invalid(double val,const char *p) {
	if(val<0) {
		fprintf(stderr,"Value of %s either invalid or not set\n",p);
		exit(1);
	}
}
