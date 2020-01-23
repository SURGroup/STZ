#include "bd_sim.hh"

double bd_sim::velocity(double x,double y) {
	double x0,y0,qu,qv;
	ls.bicubic(x,y,x0,y0);
	bicubic_velocity(x,y,qu,qv);
	return -qu*x0-qv*y0;
//	return -(fabs(x)>wallx?(x>0?wallu:-wallu):qu)*x0-qv*y0;
}

void bd_sim::bicubic_velocity(double x,double y,double &qu,double &qv) {
	const double sx=0.5,sy=0.5,s2=0.25;
	int i,j;
	double ulbx,urbx,ultx,urtx,ulby,urby,ulty,urty;
	double ulbxy,urbxy,ultxy,urtxy;
	double vlbx,vrbx,vltx,vrtx,vlby,vrby,vlty,vrty;
	double vlbxy,vrbxy,vltxy,vrtxy;
	double fx=(x-ax)*xsp,fy=(y-ay)*ysp,gx,gy;
	i=int(fx);j=int(fy);

	if(j<1) {
		gy=1-fy;
		if(i<1) {
			gx=1-fx;
			urbx=(fm[2].u-fm->u)*sx;
			urtx=(fm[me+2].u-fm[me].u)*sx;
			ulty=(fm[2*me].u-fm->u)*sy;
			urty=(fm[2*me+1].u-fm[1].u)*sy;
			urtxy=(fm[2*me+2].u-fm[2].u-fm[2*me].u+fm->u)*s2;
			qu=gy*gy*(fm->u*gx*gx+fm[1].u*fx*(2-fx)-urbx*fx*gx)
				+fy*(2-fy)*(fm[me].u*gx*gx+fm[me+1].u*fx*(2-fx)-urtx*fx*gx)
				-fy*gy*(ulty*gx*gx+urty*fx*(2-fx)-urtxy*fx*gx);
			vrbx=(fm[2].v-fm->v)*sx;
			vrtx=(fm[me+2].v-fm[me].v)*sx;
			vlty=(fm[2*me].v-fm->v)*sy;
			vrty=(fm[2*me+1].v-fm[1].v)*sy;
			vrtxy=(fm[2*me+2].v-fm[2].v-fm[2*me].v+fm->v)*s2;
			qv=gy*gy*(fm->v*gx*gx+fm[1].v*fx*(2-fx)-vrbx*fx*gx)
				+fy*(2-fy)*(fm[me].v*gx*gx+fm[me+1].v*fx*(2-fx)-vrtx*fx*gx)
				-fy*gy*(vlty*gx*gx+vrty*fx*(2-fx)-vrtxy*fx*gx);
		} else if(i>=me-2) {
			fx-=double(me-2);
			gx=1-fx;
			c_field *f=fm+me-2;
			ulbx=(f[1].u-f[-1].u)*sx;
			ultx=(f[me+1].u-f[me-1].u)*sx;
			ulty=(f[2*me].u-f->u)*sy;
			urty=(f[2*me+1].u-f[1].u)*sy;
			ultxy=(f[2*me+1].u-f[1].u-f[2*me-1].u+f[-1].u)*s2;
			qu=gy*gy*(ulbx*fx*gx+f->u*(1-fx*fx)+f[1].u*fx*fx)
				+fy*(2-fy)*(ultx*fx*gx+f[me].u*(1-fx*fx)+f[me+1].u*fx*fx)
				-fy*gy*(ultxy*fx*gx+ulty*(1-fx*fx)+urty*fx*fx);
			vlbx=(f[1].v-f[-1].v)*sx;
			vltx=(f[me+1].v-f[me-1].v)*sx;
			vlty=(f[2*me].v-f->v)*sy;
			vrty=(f[2*me+1].v-f[1].v)*sy;
			vltxy=(f[2*me+1].v-f[1].v-f[2*me-1].v+f[-1].v)*s2;
			qv=gy*gy*(vlbx*fx*gx+f->v*(1-fx*fx)+f[1].v*fx*fx)
				+fy*(2-fy)*(vltx*fx*gx+f[me].v*(1-fx*fx)+f[me+1].v*fx*fx)
				-fy*gy*(vltxy*fx*gx+vlty*(1-fx*fx)+vrty*fx*fx);
		} else {
			fx-=double(i);
			gx=1-fx;
			c_field *f=fm+i;
			ulbx=(f[1].u-f[-1].u)*sx;
			ultx=(f[me+1].u-f[me-1].u)*sx;
			urbx=(f[2].u-f->u)*sx;
			urtx=(f[me+2].u-f[me].u)*sx;
			ulty=(f[2*me].u-f->u)*sy;
			urty=(f[2*me+1].u-f[1].u)*sy;
			ultxy=(f[2*me+1].u-f[1].u-f[2*me-1].u+f[-1].u)*s2;
			urtxy=(f[2*me+2].u-f[2].u-f[2*me].u+f->u)*s2;
			qu=gy*gy*(ulbx*fx*gx*gx+f->u*gx*gx*(1+2*fx)+f[1].u*fx*fx*(3-2*fx)-urbx*fx*fx*gx)
				+fy*(2-fy)*(ultx*fx*gx*gx+f[me].u*gx*gx*(1+2*fx)+f[me+1].u*fx*fx*(3-2*fx)-urtx*fx*fx*gx)
				-fy*gy*(ultxy*fx*gx*gx+ulty*gx*gx*(1+2*fx)+urty*fx*fx*(3-2*fx)-urtxy*fx*fx*gx);
			vlbx=(f[1].v-f[-1].v)*sx;
			vltx=(f[me+1].v-f[me-1].v)*sx;
			vrbx=(f[2].v-f->v)*sx;
			vrtx=(f[me+2].v-f[me].v)*sx;
			vlty=(f[2*me].v-f->v)*sy;
			vrty=(f[2*me+1].v-f[1].v)*sy;
			vltxy=(f[2*me+1].v-f[1].v-f[2*me-1].v+f[-1].v)*s2;
			vrtxy=(f[2*me+2].v-f[2].v-f[2*me].v+f->v)*s2;
			qv=gy*gy*(vlbx*fx*gx*gx+f->v*gx*gx*(1+2*fx)+f[1].v*fx*fx*(3-2*fx)-vrbx*fx*fx*gx)
				+fy*(2-fy)*(vltx*fx*gx*gx+f[me].v*gx*gx*(1+2*fx)+f[me+1].v*fx*fx*(3-2*fx)-vrtx*fx*fx*gx)
				-fy*gy*(vltxy*fx*gx*gx+vlty*gx*gx*(1+2*fx)+vrty*fx*fx*(3-2*fx)-vrtxy*fx*fx*gx);
		}
	} else if(j>=ne-2) {
		fy-=double(ne-2);
		gy=1-fy;
		if(i<1) {
			gx=1-fx;
			c_field *f=fm+mne-2*me;
			urbx=(f[2].u-f->u)*sx;
			urtx=(f[me+2].u-f[me].u)*sx;
			ulby=(f[me].u-f[-me].u)*sy;
			urby=(f[me+1].u-f[-me+1].u)*sy;
			urbxy=(f[me+2].u-f[-me+2].u-f[me].u+f[-me].u)*s2;
			qu=fy*gy*(ulby*gx*gx+urby*fx*(2-fx)-urbxy*fx*gx)
				+(1-fy*fy)*(f->u*gx*gx+f[1].u*fx*(2-fx)-urbx*fx*gx)
				+fy*fy*(f[me].u*gx*gx+f[me+1].u*fx*(2-fx)-urtx*fx*gx);
			vrbx=(f[2].v-f->v)*sx;
			vrtx=(f[me+2].v-f[me].v)*sx;
			vlby=(f[me].v-f[-me].v)*sy;
			vrby=(f[me+1].v-f[-me+1].v)*sy;
			vrbxy=(f[me+2].v-f[-me+2].v-f[me].v+f[-me].v)*s2;
			qv=fy*gy*(vlby*gx*gx+vrby*fx*(2-fx)-vrbxy*fx*gx)
				+(1-fy*fy)*(f->v*gx*gx+f[1].v*fx*(2-fx)-vrbx*fx*gx)
				+fy*fy*(f[me].v*gx*gx+f[me+1].v*fx*(2-fx)-vrtx*fx*gx);
		} else if(i>=me-2) {
			fx-=double(me-2);
			gx=1-fx;
			c_field *f=fm+mne-me-2;
			ulbx=(f[1].u-f[-1].u)*sx;
			ultx=(f[me+1].u-f[me-1].u)*sx;
			ulby=(f[me].u-f[-me].u)*sy;
			urby=(f[me+1].u-f[-me+1].u)*sy;
			ulbxy=(f[me+1].u-f[-me+1].u-f[me-1].u+f[-me-1].u)*s2;
			qu=fy*gy*(ulbxy*fx*gx+ulby*(1-fx*fx)+urby*fx*fx)
				+(1-fy*fy)*(ulbx*fx*gx+f->u*(1-fx*fx)+f[1].u*fx*fx)
				+fy*fy*(ultx*fx*gx+f[me].u*(1-fx*fx)+f[me+1].u*fx*fx);
			vlbx=(f[1].v-f[-1].v)*sx;
			vltx=(f[me+1].v-f[me-1].v)*sx;
			vlby=(f[me].v-f[-me].v)*sy;
			vrby=(f[me+1].v-f[-me+1].v)*sy;
			vlbxy=(f[me+1].v-f[-me+1].v-f[me-1].v+f[-me-1].v)*s2;
			qv=fy*gy*(vlbxy*fx*gx+vlby*(1-fx*fx)+vrby*fx*fx)
				+(1-fy*fy)*(vlbx*fx*gx+f->v*(1-fx*fx)+f[1].v*fx*fx)
				+fy*fy*(vltx*fx*gx+f[me].v*(1-fx*fx)+f[me+1].v*fx*fx);
		} else {
			fx-=double(i);
			gx=1-fx;
			c_field *f=fm+mne-2*me+i;
			ulbx=(f[1].u-f[-1].u)*sx;
			ultx=(f[me+1].u-f[me-1].u)*sx;
			urbx=(f[2].u-f->u)*sx;
			urtx=(f[me+2].u-f[me].u)*sx;
			ulby=(f[me].u-f[-me].u)*sy;
			urby=(f[me+1].u-f[-me+1].u)*sy;
			ulbxy=(f[me+1].u-f[-me+1].u-f[me-1].u+f[-me-1].u)*s2;
			urbxy=(f[me+2].u-f[-me+2].u-f[me].u+f[-me].u)*s2;
			qu=fy*gy*(ulbxy*fx*gx*gx+ulby*gx*gx*(1+2*fx)+urby*fx*fx*(3-2*fx)-urbxy*fx*fx*gx)
				+(1-fy*fy)*(ulbx*fx*gx*gx+f->u*gx*gx*(1+2*fx)+f[1].u*fx*fx*(3-2*fx)-urbx*fx*fx*gx)
				+fy*fy*(ultx*fx*gx*gx+f[me].u*gx*gx*(1+2*fx)+f[me+1].u*fx*fx*(3-2*fx)-urtx*fx*fx*gx);
			vlbx=(f[1].v-f[-1].v)*sx;
			vltx=(f[me+1].v-f[me-1].v)*sx;
			vrbx=(f[2].v-f->v)*sx;
			vrtx=(f[me+2].v-f[me].v)*sx;
			vlby=(f[me].v-f[-me].v)*sy;
			vrby=(f[me+1].v-f[-me+1].v)*sy;
			vlbxy=(f[me+1].v-f[-me+1].v-f[me-1].v+f[-me-1].v)*s2;
			vrbxy=(f[me+2].v-f[-me+2].v-f[me].v+f[-me].v)*s2;
			qv=fy*gy*(vlbxy*fx*gx*gx+vlby*gx*gx*(1+2*fx)+vrby*fx*fx*(3-2*fx)-vrbxy*fx*fx*gx)
				+(1-fy*fy)*(vlbx*fx*gx*gx+f->v*gx*gx*(1+2*fx)+f[1].v*fx*fx*(3-2*fx)-vrbx*fx*fx*gx)
				+fy*fy*(vltx*fx*gx*gx+f[me].v*gx*gx*(1+2*fx)+f[me+1].v*fx*fx*(3-2*fx)-vrtx*fx*fx*gx);
		}
	} else {
		fy-=double(j);
		gy=1-fy;
		if(i<1) {
			gx=1-fx;
			c_field *f=fm+j*me;
			urbx=(f[2].u-f->u)*sx;
			urtx=(f[me+2].u-f[me].u)*sx;
			ulby=(f[me].u-f[-me].u)*sy;
			urby=(f[me+1].u-f[-me+1].u)*sy;
			ulty=(f[2*me].u-f->u)*sy;
			urty=(f[2*me+1].u-f[1].u)*sy;
			urbxy=(f[me+2].u-f[-me+2].u-f[me].u+f[-me].u)*s2;
			urtxy=(f[2*me+2].u-f[2].u-f[2*me].u+f->u)*s2;
			qu=fy*gy*gy*(ulby*gx*gx+urby*fx*(2-fx)-urbxy*fx*gx)
				+gy*gy*(1+2*fy)*(f->u*gx*gx+f[1].u*fx*(2-fx)-urbx*fx*gx)
				+fy*fy*(3-2*fy)*(f[me].u*gx*gx+f[me+1].u*fx*(2-fx)-urtx*fx*gx)
				-fy*fy*gy*(ulty*gx*gx+urty*fx*(2-fx)-urtxy*fx*gx);
			vrbx=(f[2].v-f->v)*sx;
			vrtx=(f[me+2].v-f[me].v)*sx;
			vlby=(f[me].v-f[-me].v)*sy;
			vrby=(f[me+1].v-f[-me+1].v)*sy;
			vlty=(f[2*me].v-f->v)*sy;
			vrty=(f[2*me+1].v-f[1].v)*sy;
			vrbxy=(f[me+2].v-f[-me+2].v-f[me].v+f[-me].v)*s2;
			vrtxy=(f[2*me+2].v-f[2].v-f[2*me].v+f->v)*s2;
			qv=fy*gy*gy*(vlby*gx*gx+vrby*fx*(2-fx)-vrbxy*fx*gx)
				+gy*gy*(1+2*fy)*(f->v*gx*gx+f[1].v*fx*(2-fx)-vrbx*fx*gx)
				+fy*fy*(3-2*fy)*(f[me].v*gx*gx+f[me+1].v*fx*(2-fx)-vrtx*fx*gx)
				-fy*fy*gy*(vlty*gx*gx+vrty*fx*(2-fx)-vrtxy*fx*gx);
		} else if(i>=me-2) {
			fx-=double(me-2);
			gx=1-fx;
			c_field *f=fm+(j+1)*me-2;
			ulbx=(f[1].u-f[-1].u)*sx;
			ultx=(f[me+1].u-f[me-1].u)*sx;
			ulby=(f[me].u-f[-me].u)*sy;
			urby=(f[me+1].u-f[-me+1].u)*sy;
			ulty=(f[2*me].u-f->u)*sy;
			urty=(f[2*me+1].u-f[1].u)*sy;
			ulbxy=(f[me+1].u-f[-me+1].u-f[me-1].u+f[-me-1].u)*s2;
			ultxy=(f[2*me+1].u-f[1].u-f[2*me-1].u+f[-1].u)*s2;
			qu=fy*gy*gy*(ulbxy*fx*gx+ulby*(1-fx*fx)+urby*fx*fx)
				+gy*gy*(1+2*fy)*(ulbx*fx*gx+f->u*(1-fx*fx)+f[1].u*fx*fx)
				+fy*fy*(3-2*fy)*(ultx*fx*gx+f[me].u*(1-fx*fx)+f[me+1].u*fx*fx)
				-fy*fy*gy*(ultxy*fx*gx+ulty*(1-fx*fx)+urty*fx*fx);
			vlbx=(f[1].v-f[-1].v)*sx;
			vltx=(f[me+1].v-f[me-1].v)*sx;
			vlby=(f[me].v-f[-me].v)*sy;
			vrby=(f[me+1].v-f[-me+1].v)*sy;
			vlty=(f[2*me].v-f->v)*sy;
			vrty=(f[2*me+1].v-f[1].v)*sy;
			vlbxy=(f[me+1].v-f[-me+1].v-f[me-1].v+f[-me-1].v)*s2;
			vltxy=(f[2*me+1].v-f[1].v-f[2*me-1].v+f[-1].v)*s2;
			qv=fy*gy*gy*(vlbxy*fx*gx+vlby*(1-fx*fx)+vrby*fx*fx)
				+gy*gy*(1+2*fy)*(vlbx*fx*gx+f->v*(1-fx*fx)+f[1].v*fx*fx)
				+fy*fy*(3-2*fy)*(vltx*fx*gx+f[me].v*(1-fx*fx)+f[me+1].v*fx*fx)
				-fy*fy*gy*(vltxy*fx*gx+vlty*(1-fx*fx)+vrty*fx*fx);
		} else {
			fx-=double(i);
			gx=1-fx;
			c_field *f=fm+j*me+i;
			ulbx=(f[1].u-f[-1].u)*sx;
			ultx=(f[me+1].u-f[me-1].u)*sx;
			urbx=(f[2].u-f->u)*sx;
			urtx=(f[me+2].u-f[me].u)*sx;
			ulby=(f[me].u-f[-me].u)*sy;
			urby=(f[me+1].u-f[-me+1].u)*sy;
			ulty=(f[2*me].u-f->u)*sy;
			urty=(f[2*me+1].u-f[1].u)*sy;
			ulbxy=(f[me+1].u-f[-me+1].u-f[me-1].u+f[-me-1].u)*s2;
			ultxy=(f[2*me+1].u-f[1].u-f[2*me-1].u+f[-1].u)*s2;
			urbxy=(f[me+2].u-f[-me+2].u-f[me].u+f[-me].u)*s2;
			urtxy=(f[2*me+2].u-f[2].u-f[2*me].u+f->u)*s2;
			qu=fy*gy*gy*(ulbxy*fx*gx*gx+ulby*gx*gx*(1+2*fx)+urby*fx*fx*(3-2*fx)-urbxy*fx*fx*gx)
				+gy*gy*(1+2*fy)*(ulbx*fx*gx*gx+f->u*gx*gx*(1+2*fx)+f[1].u*fx*fx*(3-2*fx)-urbx*fx*fx*gx)
				+fy*fy*(3-2*fy)*(ultx*fx*gx*gx+f[me].u*gx*gx*(1+2*fx)+f[me+1].u*fx*fx*(3-2*fx)-urtx*fx*fx*gx)
				-fy*fy*gy*(ultxy*fx*gx*gx+ulty*gx*gx*(1+2*fx)+urty*fx*fx*(3-2*fx)-urtxy*fx*fx*gx);
			vlbx=(f[1].v-f[-1].v)*sx;
			vltx=(f[me+1].v-f[me-1].v)*sx;
			vrbx=(f[2].v-f->v)*sx;
			vrtx=(f[me+2].v-f[me].v)*sx;
			vlby=(f[me].v-f[-me].v)*sy;
			vrby=(f[me+1].v-f[-me+1].v)*sy;
			vlty=(f[2*me].v-f->v)*sy;
			vrty=(f[2*me+1].v-f[1].v)*sy;
			vlbxy=(f[me+1].v-f[-me+1].v-f[me-1].v+f[-me-1].v)*s2;
			vltxy=(f[2*me+1].v-f[1].v-f[2*me-1].v+f[-1].v)*s2;
			vrbxy=(f[me+2].v-f[-me+2].v-f[me].v+f[-me].v)*s2;
			vrtxy=(f[2*me+2].v-f[2].v-f[2*me].v+f->v)*s2;
			qv=fy*gy*gy*(vlbxy*fx*gx*gx+vlby*gx*gx*(1+2*fx)+vrby*fx*fx*(3-2*fx)-vrbxy*fx*fx*gx)
				+gy*gy*(1+2*fy)*(vlbx*fx*gx*gx+f->v*gx*gx*(1+2*fx)+f[1].v*fx*fx*(3-2*fx)-vrbx*fx*fx*gx)
				+fy*fy*(3-2*fy)*(vltx*fx*gx*gx+f[me].v*gx*gx*(1+2*fx)+f[me+1].v*fx*fx*(3-2*fx)-vrtx*fx*fx*gx)
				-fy*fy*gy*(vltxy*fx*gx*gx+vlty*gx*gx*(1+2*fx)+vrty*fx*fx*(3-2*fx)-vrtxy*fx*fx*gx);
		}
	}
}
