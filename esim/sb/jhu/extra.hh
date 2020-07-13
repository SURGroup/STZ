#ifndef EXTRA_HH
#define EXTRA_HH

#include "shear_sim.hh"

/** Checks to see if two strings are equal.
 * \param[in] (p1,p2) the two strings.
 * \return True if the are equal, false if they are not. */
inline bool se(const char* p1,const char *p2) {
    return strcmp(p1,p2)==0;
}

void read_chi_from_file(const char *filename,shear_sim &sim,double u0,double uscale,double TZ);
bool read_arg(char **argv,const char *p2,int &i,int argc,double &val,const char *name,const char *unit);

#endif
