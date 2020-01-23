#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "gp_matrix.hh"

int main(int argc,char **argv) {

	// Check for the correct number of command-line arguments
	if(argc!=5) {
		fputs("Usage: ./y_transect <output_dir> <suffix> <x_pos> <n_max>\n",stderr);
		return 1;
	}

	// Read in all of the frames
	char buf[256];
	int i,j,nmax=atoi(argv[4]),m=0,n=0;
	float xs=atof(argv[3]),**yt=new float*[nmax+1],ay=0,dy=0;
	for(i=0;i<=nmax;i++) {

		// Assemble the filename of the frame and read it in
		sprintf(buf,"%s/%s.%d",argv[1],argv[2],i);
		gp_matrix gp(buf);

		// On the first frame, store the dimension. On subsequent frames, check
		// that the dimensions agree
		if(i==0) {
			m=gp.m;n=gp.n;ay=*(gp.y);dy=gp.dx;
		} else if(m!=gp.m||n!=gp.n) {
			fprintf(stderr,"Size mismatch in frame %d\n",i);
			return 1;
		}

		// Allocate memory and store the transcript
		yt[i]=new float[n];
		gp.y_transect(xs,yt[i]);
	}

	// Save the output file
	sprintf(buf,"%s_%s.yt_%s",argv[1],argv[2],argv[3]);
	FILE *fp=fopen(buf,"w");
	if(fp==NULL) {
		fputs("Can't open output file\n",stderr);
		return 1;
	}
	for(j=0;j<n;j++) {
		fprintf(fp,"%g",ay+j*dy);
		for(i=0;i<=nmax;i++) fprintf(fp," %g",yt[i][j]);
		fputs("\n",fp);
	}
	fclose(fp);

	// Deallocate memory
	for(i=0;i<=nmax;i++) delete [] yt[i];
	delete [] yt;
}
