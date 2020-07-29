#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc<5||argc>6) {
		fputs("Usage: ./colorbar <output_png> <width> <height> <palette> [gamma]\n",stderr);
		return 1;
	}

	// Check the dimensions are sensible
	int wid=atoi(argv[2]),hei=atoi(argv[3]),ptype=0;
	if(wid<=0||wid>16777216||hei<=0||hei>16777216) {
		fputs("Image dimensions out of bounds\n",stderr);
		return 1;
	}

	// Check the palette number
    ptype=atoi(argv[4]);
    if(ptype<0||ptype>7) {
        fputs("Palette number out of range\n",stderr);
        return 1;
    }

    // Read the gamma value if specified
    double gamma=1.;
    if(argc==6) {
        gamma=atof(argv[5]);
        if(gamma<=0||gamma>1e6) {
            fputs("Gamma value it outside a typical range\n",stderr);
            return 1;
        }
    }

	// Create and output the color bar
	gp_matrix_bitmap_colorbar(argv[1],wid,hei,ptype,gamma);
}
