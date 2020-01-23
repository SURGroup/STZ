#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int nframes=500;
float *phi,*xp,*yp,*f;

int main(int argc,char **argv) {
	if(argc!=2) {
		cerr << "Usage: ./wallforce <directory>" << endl;
		return 1;
	}

	int m,n,an=0,amn,i,j,frame,wframe,ds=sizeof(float),nwi,pwi;
	float ifloat,nz,pz,wp,zz,bnb,bnt,bpb,bpt,fd,fu,phid,phiu,nwf,pwf;
	char *floatd;
	static char fname[256];
	floatd=(char*)&ifloat;

	ifstream inf,inf2,wallf;
	ofstream wff;
	sprintf(fname,"%s/wall",argv[1]);
	wallf.open(fname,fstream::in);
	sprintf(fname,"%s.wf",argv[1]);
	wff.open(fname,fstream::out|fstream::trunc);

	for(frame=0;frame<=nframes;frame++) {
		cout << frame << endl;

		wallf >> wframe >> wp;
		if(wframe!=frame) {
			cerr << "Wall frame mismatch" << endl;
			return 1;
		}

		sprintf(fname,"%s/phi.%d",argv[1],frame);
		inf.open(fname,fstream::in|fstream::binary);

		inf.read(floatd,ds);
		m=int(ifloat+0.5);

		// Allocate memory for phi and grid
		if(frame==0) {
			an=1+(m>>1);
			amn=m*an;
			phi=new float[amn];
			xp=new float[m];
			yp=new float[an];
			f=new float[amn];
		}

		nwi=0;pwi=0;
		for(i=0;i<m;i++) {
			inf.read(floatd,ds);
			if(ifloat<-wp) nwi=i;
			if(ifloat<wp) pwi=i;
			xp[i]=ifloat;
		}

		if(pwi==m-1||nwi==m-1) {
			cerr << "Couldn't find wall intersection point" << endl;return 1;
		}

		nz=(xp[nwi]+wp)/(xp[nwi]-xp[nwi+1]);
		pz=(xp[pwi]-wp)/(xp[pwi]-xp[pwi+1]);

		n=0;
		inf.read(floatd,ds);
		while(!inf.eof()) {
			yp[n]=ifloat;
			for(i=0;i<m;i++) {
				inf.read(floatd,ds);
				phi[n*m+i]=ifloat;
			}
			n++;
			if(n>an) {
				cerr << "Not enough vertical memory" << endl;return 1;
			}
			inf.read(floatd,ds);
		}
		inf.close();

		sprintf(fname,"%s/p.%d",argv[1],frame);
		inf.open(fname,fstream::in|fstream::binary);
		sprintf(fname,"%s/s.%d",argv[1],frame);
		inf2.open(fname,fstream::in|fstream::binary);

		for(i=0;i<m+1;i++) inf.read(floatd,ds);
		for(i=0;i<m+1;i++) inf2.read(floatd,ds);
		for(j=0;j<n;j++) {
			inf.read(floatd,ds);
			inf2.read(floatd,ds);
			for(i=0;i<m;i++) {
				inf.read(floatd,ds);
				f[j*m+i]=ifloat;
				inf2.read(floatd,ds);
				f[j*m+i]-=ifloat;
			}
		}
		inf.close();
		inf2.close();

		bnb=bnt=bpb=bpt=-1;

		nwf=pwf=0;
		for(j=0;j<n-1;j++) {
			f[m*j+nwi]=2*f[m*j+nwi+1]-f[m*j+nwi+2];
			phid=(1-nz)*phi[j*m+nwi]+nz*phi[j*m+nwi+1];
			phiu=(1-nz)*phi[(j+1)*m+nwi]+nz*phi[(j+1)*m+nwi+1];
			if(phid<0) {
				fd=f[m*j+nwi]*(1-nz)+f[m*j+nwi+1]*nz;
				fu=f[m*(j+1)+nwi]*(1-nz)+f[m*(j+1)+nwi+1]*nz;
				if(phiu<0) {
					nwf+=0.5*(fd+fu)*(yp[j+1]-yp[j]);
				} else {
					zz=phid/(phid-phiu);
					fu=fd*(1-zz)+fu*zz;
					bnt=yp[j]*(1-zz)+yp[j+1]*zz;
					nwf+=0.5*(fd+fu)*(bnt-yp[j]);
				}
			} else {
				if(phiu<0) {
					fd=f[m*j+nwi]*(1-nz)+f[m*j+nwi+1]*nz;
					fu=f[m*(j+1)+nwi]*(1-nz)+f[m*(j+1)+nwi+1]*nz;
					zz=phid/(phid-phiu);
					fd=fd*(1-zz)+fu*zz;
					bnb=yp[j]*(1-zz)+yp[j+1]*zz;
					nwf+=0.5*(fd+fu)*(yp[j+1]-bnb);
				}
			}
			f[m*j+pwi+1]=2*f[m*j+pwi]-f[m*j+pwi-1];
			phid=(1-pz)*phi[j*m+pwi]+pz*phi[j*m+pwi+1];
			phiu=(1-pz)*phi[(j+1)*m+pwi]+pz*phi[(j+1)*m+pwi+1];
			if(phid<0) {
				fd=f[m*j+pwi]*(1-pz)+f[m*j+pwi+1]*pz;
				fu=f[m*(j+1)+pwi]*(1-pz)+f[m*(j+1)+pwi+1]*pz;
				if(phiu<0) {
					pwf+=0.5*(fd+fu)*(yp[j+1]-yp[j]);
				} else {
					zz=phid/(phid-phiu);
					fu=fd*(1-zz)+fu*zz;
					bpt=yp[j]*(1-zz)+yp[j+1]*zz;
					pwf+=0.5*(fd+fu)*(bpt-yp[j]);
				}
			} else {
				if(phiu<0) {
					fd=f[m*j+pwi]*(1-pz)+f[m*j+pwi+1]*pz;
					fu=f[m*(j+1)+pwi]*(1-pz)+f[m*(j+1)+pwi+1]*pz;
					zz=phid/(phid-phiu);
					fd=fd*(1-zz)+fu*zz;
					bpb=yp[j]*(1-zz)+yp[j+1]*zz;
					pwf+=0.5*(fd+fu)*(yp[j+1]-bpb);
				}
			}
		}

		wff << frame << " " << bnb << " " << bnt << " " << nwf << " " << bpb << " " << bpt << " " << pwf << endl;

	}

	delete [] phi;
	delete [] xp;
	delete [] yp;
	delete [] f;
}
