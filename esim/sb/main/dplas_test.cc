#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "stz_model.hh"
#include "fileinfo.hh"

int main(int argc,char **argv) {
	if(argc!=2) {
		fputs("Syntax: ./dplas_test <input_file>\n",stderr);
		return 1;
	}
	fileinfo fin(argv[argc-1]);
	double x,y,dc1,dc2;
	for(x=0.8;x<2;x+=0.001) {
		printf("%g",x);
		for(y=0.25;y<0.46;y+=0.05) printf(" %g",fin.stz->Dplastic(x,y,dc1,dc2));
		puts("");
	}
}
