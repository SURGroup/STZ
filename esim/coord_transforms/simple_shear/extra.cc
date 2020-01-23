#include "extra.hh"
#include "common.hh"

#include <cstdio>
#include <cstdlib>

/** Sets the chi field in a shear_sim object based on potential energy values stored in a file.
 * \param[in] filename the file to read from.
 * \param[in,out] sim a reference to the shear_sim object.
 * \param[in] u0 the displacement to apply to the energy values.
 * \param[in] uscale the scale to apply to the energy values.
 * \param[in] TZ the STZ formation energy. */
void read_chi_from_file(const char *filename,shear_sim &sim,double u0,double uscale,double TZ) {
    FILE *fp=safe_fopen(filename,"r");
    int mm,nn;
    double chi_base,chi_scale;
    if(fscanf(fp,"%d %d %lg %lg",&mm,&nn,&chi_base,&chi_scale)!=4) {
        fputs("Error reading header information\n",stderr);
        exit(1);
    }

    // Check that the dimensions make sense, and allocate memory
    if(mm<=0 || nn<=0) {
        fputs("Grid dimensions are invalid\n",stderr);
        exit(1);
    }
    double *f=new double[mm*nn];

    // Read in the energy values from the file
    for(int j=0;j<nn;j++) for(int i=0;i<mm;i++) {
        if(fscanf(fp,"%lg",f+i+mm*j)!=1) {
            fputs("Error reading energy information\n",stderr);
            exit(1);
        }
    }

    // Map energy values onto effective temperature, chi
    for(int j=0;j<nn;j++) for(int i=0;i<mm;i++)
        f[i+mm*j]=(f[i+mm*j]-u0)*uscale;

    printf("# chi grid dimensions : %d by %d\n"
           "# chi base scale      : %g K\n"
           "# chi range scale     : %g K\n",mm,nn,chi_base,chi_scale);
    sim.initialize_chi_bicubic(mm,nn,f,chi_base/TZ,chi_scale/TZ);

    // Free the dynamically allocated memory
    delete [] f;
}

/** Reads in a command-line argument and sets a corresponding parameter.
 * \param[in] argv a pointer to the command-line arguments.
 * \param[in] p2 the string describing the parameter.
 * \param[in] i the index of the command-line argument to look at.
 * \param[in,out] val the value to set.
 * \param[in] name the longer name of the parameter.
 * \param[in] unit the unit of the parameter.
 * \return True if a parameter was set, false otherwise. */
bool read_arg(char **argv,const char *p2,int &i,int argc,double &val,const char *name,const char *unit) {
    if(se(argv[i],p2)) {
        if(++i==argc) {
            fputs("Error reading command-line arguments\n",stderr);
            exit(1);
        }
        printf("The %s is set to %g%s\n",name,val=atof(argv[i]),unit);
        return true;
    }
    return false;
}
