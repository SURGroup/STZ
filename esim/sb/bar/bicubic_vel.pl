#!/usr/bin/perl
open B,">bicubic_vel.cc";

# Print header
print B <<EOF;
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

EOF

# Table of bicubic interpolation function components
@func=(["","gx*gx","fx*(2-fx)","-fx*gx"],
	["fx*gx","(1-fx*fx)","fx*fx",""],
	["fx*gx*gx","gx*gx*(1+2*fx)","fx*fx*(3-2*fx)","-fx*fx*gx"]);

# Table of bicubic interpolation coefficients
@expr=(["plbxy","plby","prby","prbxy"],["plbx","f->p","f[1].p","prbx"],
	["pltx","f[me].p","f[me+1].p","prtx"],["pltxy","plty","prty","prtxy"]);

# Assemble strings
foreach $k (0..2) {
	foreach $l (0..3) {
		$_=$func[$k][$l];
		$sgn[$k][$l]=m/^-/?"-":"+";
		$nsgn[$k][$l]=m/^-/?"-":"";
		s/^-//g;$e[$k][$l]=$_;
		s/x/y/g;$f[$k][$l]=$_;
	}
}

foreach $b (0,1,2) {

	if ($b==0) {
		print B "\tif(j<1) {\n";
	} elsif ($b==1) {
		print B "\t} else if(j>=ne-2) {\n";
		print B "\t\tfy-=double(ne-2);\n";
	} else {
		print B "\t} else {\n";
		print B "\t\tfy-=double(j);\n";
	}
	print B "\t\tgy=1-fy;\n";

	foreach $a (0,1,2) {

		if ($a==0) {
			print B "\t\tif(i<1) {\n";
		} elsif($a==1) {
			print B "\t\t} else if(i>=me-2) {\n";
			print B "\t\t\tfx-=double(me-2);\n";
		} else {
			print B "\t\t} else {\n";
			print B "\t\t\tfx-=double(i);\n";
		}
		print B "\t\t\tgx=1-fx;\n";

		$c=3*$b+$a;
		if ($c==1) {$ij="me-2";}
		elsif ($c==2) {$ij="i";}
		elsif ($c==3) {$ij="mne-2*me";}
		elsif ($c==4) {$ij="mne-me-2";}
		elsif ($c==5) {$ij="mne-2*me+i";}
		elsif ($c==6) {$ij="j*me";}
		elsif ($c==7) {$ij="(j+1)*me-2";}
		elsif ($c==8) {$ij="j*me+i";}

		print B "\t\t\tc_field *f=fm+$ij;\n" unless $c==0;

		$_="";
		if($a!=0) {
			$_.="\t\t\tplbx=(f[1].p-f[-1].p)*sx;\n";
			$_.="\t\t\tpltx=(f[me+1].p-f[me-1].p)*sx;\n";
		}
		if($a!=1) {
			$_.="\t\t\tprbx=(f[2].p-f->p)*sx;\n";
			$_.="\t\t\tprtx=(f[me+2].p-f[me].p)*sx;\n";
		}
		if($b!=0) {
			$_.="\t\t\tplby=(f[me].p-f[-me].p)*sy;\n";
			$_.="\t\t\tprby=(f[me+1].p-f[-me+1].p)*sy;\n";
		}
		if($b!=1) {
			$_.="\t\t\tplty=(f[2*me].p-f->p)*sy;\n";
			$_.="\t\t\tprty=(f[2*me+1].p-f[1].p)*sy;\n";
		}
		if($b!=0&&$a!=0) {
			$_.="\t\t\tplbxy=(f[me+1].p-f[-me+1].p-f[me-1].p+f[-me-1].p)*s2;\n";
		}
		if($b!=1&&$a!=0) {
			$_.="\t\t\tpltxy=(f[2*me+1].p-f[1].p-f[2*me-1].p+f[-1].p)*s2;\n";
		}
		if($b!=0&&$a!=1) {
			$_.="\t\t\tprbxy=(f[me+2].p-f[-me+2].p-f[me].p+f[-me].p)*s2;\n";
		}
		if($b!=1&&$a!=1) {
			$_.="\t\t\tprtxy=(f[2*me+2].p-f[2].p-f[2*me].p+f->p)*s2;\n";
		}
		s/f/fm/g if $c==0;

		$st=0;
		foreach $l (0..3) {
			next if $f[$b][$l] eq "";

			if ($st==0) {
				$_.="\t\t\tqp=$nsgn[$b][$l]";$st=1;
			} else {
				$_.="\n\t\t\t\t$sgn[$b][$l]";
			}
			$_.="$f[$b][$l]*(";

			$st2=0;
			foreach $k (0..3) {
				next if $e[$a][$k] eq "";

				if($st2==0) {
					$_.="$nsgn[$a][$k]";$st2=1;
				} else {
					$_.="$sgn[$a][$k]";
				}

				$_.="$expr[$l][$k]*$e[$a][$k]";
			}
			$_.=")";

		}
		$_.=";\n";
		s/f->/fm->/g, s/f\[/fm\[/g if $c==0;

		s/p/u/g;print B;
		s/u/v/g;print B;
	}
	print B "\t\t}\n";
}
print B "\t}\n}\n";
