#include "common.hh"
#include "common.hh"

FILE* safe_fopen(const char* filename,const char* mode) {
	FILE *temp=fopen(filename,mode);
	if(temp==NULL) fprintf(stderr,"sbsim: error opening file \"%s\"",filename);
	return temp;
}

void fatal_error(const char *p,int code) {
	fprintf(stderr,"sbsim: %s\n",p);
	exit(code);
}
