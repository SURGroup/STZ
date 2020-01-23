#include <cstdio>
#include <cstdlib>
#include "gp_matrix.hh"

int main(int argc,char **argv) {
	int i,k;

	if(argc!=4) {
		fputs("Syntax: ./x_sequence <prefix> <field> <output>\n",stderr);
		return 1;
	}

	char buf[256];
	float xti[51][1025],phit[1025];
	sprintf(buf,"%s/%s.0",argv[1],argv[2]);
	gp_matrix gp(buf);
	sprintf(buf,"%s/phi.0",argv[1]);
	gp_matrix phi(buf);
	if(gp.m!=1025||phi.m!=1025) {
		fputs("Size error\n",stderr);
		return 1;
	}
	gp.x_transect(0,*xti);
	phi.x_transect(0,phit);
	for(i=0;i<1025;i++) if(phit[i]>0) xti[0][i]=0;

	for(k=1;k<=50;k++) {
		sprintf(buf,"%s/%s.%d",argv[1],argv[2],k);
		gp.read(buf);
		sprintf(buf,"%s/phi.%d",argv[1],k);
		phi.read(buf);
		gp.x_transect(0,xti[k]);
		phi.x_transect(0,phit);
		for(i=0;i<1025;i++) if(phit[i]>0) xti[k][i]=0;
	}

	FILE *fp=fopen(argv[3],"w");
	if(fp==NULL) {
		fputs("File output error\n",stderr);
		return 1;
	}

	for(i=480;i<1025;i++) {
		fprintf(fp,"%g",gp.x[i]);
		for(k=0;k<=50;k++) fprintf(fp," %g",xti[k][i]);
		fputs("\n",fp);
	}
}
