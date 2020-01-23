open B,">bicubic2.cc";

@func=(["","mx*mx","x*(2-x)","-x*mx"],
	["x*mx","(1-x*x)","x*x",""],
	["x*mx*mx","mx*mx*(1+2*x)","x*x*(3-2*x)","-x*x*mx"]);

@dfunc=(["","-2*mx","2*mx","(2*x-1)"],
	["(1-2*x)","-2*x","2*x",""],
	["mx*(1-3*x)","-6*x*mx","6*x*mx","x*(3*x-2)"]);

@expr=(["plbxy","plby","prby","prbxy"],["plbx","fe[ij]","fe[ij+1]","prbx"],
	["pltx","fe[ij+m]","fe[ij+m+1]","prtx"],["pltxy","plty","prty","prtxy"]);

foreach $k (0..2) {
	foreach $l (0..3) {
		$_=$func[$k][$l];
		$sgn[$k][$l]=m/^-/?"-":"+";
		$nsgn[$k][$l]=m/^-/?"-":"";
		s/^-//g;
		s/x/e/g;$e[$k][$l]=$_;
		s/e/f/g;$f[$k][$l]=$_;

		$_=$dfunc[$k][$l];
		$dsgn[$k][$l]=m/^-/?"-":"+";
		$ndsgn[$k][$l]=m/^-/?"-":"";
		s/^-//g;
		s/x/e/g;$de[$k][$l]=$_;
		s/e/f/g;$df[$k][$l]=$_;
	}
}

foreach $b (0,1,2) {

	if ($b==0) {
		print B "\tif(j<1) {\n";
	} elsif ($b==1) {
		print B "\t} else if(j>=n-2) {\n";
		print B "\t\tf-=double(n-2);\n";
	} else {
		print B "\t} else {\n";
		print B "\t\tf-=double(j);\n";
	}
	print B "\t\tmf=1-f;\n";

	foreach $a (0,1,2) {

		if ($a==0) {
			print B "\t\tif(i<1) {\n";
		} elsif($a==1) {
			print B "\t\t} else if(i>=m-2) {\n";
			print B "\t\t\te-=double(m-2);\n";
		} else {
			print B "\t\t} else {\n";
			print B "\t\t\te-=double(i);\n";
		}
		print B "\t\t\tme=1-e;\n";

		$c=3*$b+$a;
		if ($c==0) {$ij="0";}
		elsif ($c==1) {$ij="m-2";}
		elsif ($c==2) {$ij="i";}
		elsif ($c==3) {$ij="mn-2*m";}
		elsif ($c==4) {$ij="mn-m-2";}
		elsif ($c==5) {$ij="mn-2*m+i";}
		elsif ($c==6) {$ij="j*m";}
		elsif ($c==7) {$ij="(j+1)*m-2";}
		elsif ($c==8) {$ij="j*m+i";}

		print B "\t\t\tij=$ij;\n";

		if($a!=0) {
			print B "\t\t\tplbx=(fe[ij+1]-fe[ij-1])*dxf;\n";
			print B "\t\t\tpltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;\n";
		}
		if($a!=1) {
			print B "\t\t\tprbx=(fe[ij+2]-fe[ij])*dxf;\n";
			print B "\t\t\tprtx=(fe[ij+m+2]-fe[ij+m])*dxf;\n";
		}
		if($b!=0) {
			print B "\t\t\tplby=(fe[ij+m]-fe[ij-m])*dyf;\n";
			print B "\t\t\tprby=(fe[ij+m+1]-fe[ij-m+1])*dyf;\n";
		}
		if($b!=1) {
			print B "\t\t\tplty=(fe[ij+2*m]-fe[ij])*dyf;\n";
			print B "\t\t\tprty=(fe[ij+2*m+1]-fe[ij+1])*dyf;\n";
		}
		if($b!=0&&$a!=0) {
			print B "\t\t\tplbxy=(fe[ij+m+1]-fe[ij-m+1]-fe[ij+m-1]+fe[ij-m-1])*d2f;\n";
		}
		if($b!=1&&$a!=0) {
			print B "\t\t\tpltxy=(fe[ij+2*m+1]-fe[ij+1]-fe[ij+2*m-1]+fe[ij-1])*d2f;\n";
		}
		if($b!=0&&$a!=1) {
			print B "\t\t\tprbxy=(fe[ij+m+2]-fe[ij-m+2]-fe[ij+m]+fe[ij-m])*d2f;\n";
		}
		if($b!=1&&$a!=1) {
			print B "\t\t\tprtxy=(fe[ij+2*m+2]-fe[ij+2]-fe[ij+2*m]+fe[ij])*d2f;\n";
		}

		$st=0;
		foreach $l (0..3) {
			next if $f[$b][$l] eq "";

			if ($st==0) {
				print B "\t\t\treturn $nsgn[$b][$l]";$st=1;
			} else {
				print B "\n\t\t\t\t$sgn[$b][$l]";
			}
			print B "$f[$b][$l]*(";

			$st2=0;
			foreach $k (0..3) {
				next if $e[$a][$k] eq "";

				if($st2==0) {
					print B "$nsgn[$a][$k]";$st2=1;
				} else {
					print B "$sgn[$a][$k]";
				}

				print B "$expr[$l][$k]*$e[$a][$k]";
			}
			print B ")";

		}
		print B ";\n";
	}
	print B "\t\t}\n";
}
print B "\t}\n";
