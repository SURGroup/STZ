#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "stz_model.hh"
#include "fileinfo.hh"

int main(int argc,char **argv) {
	if(argc!=2&&(argc!=3||strcmp(argv[1],"-h")!=0)) {
		fputs("Syntax: ./crack_digest [-h] <input_file>\n",stderr);
		return 1;
	}
	fileinfo fin(argv[argc-1]);
	if(argc==2) fin.digest();
	else fin.output_html();
}
