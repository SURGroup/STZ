open E,"frac_runs.csv" or die "Can't open table\n";
open F,">crack_table.html";

%mon=("Jan",1,"Feb",2,"Mar",3,"Apr",4,"May",5,"Jun",6,"Jul",7,"Aug",8,"Sep",9,"Oct",10,"Nov",11,"Dec",12);
@md=(31,31,28,31,30,31,30,31,31,30,31,30);
$ed[1]=0;$ed[$_]=$ed[$_-1]+$md[$_-1] foreach 2..12;

print F <<EOF;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN"
"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
	<title>STZ crack simulation table</title>
<!--	<base href="http://math.berkeley.edu/~chr/sb/" /> -->
	<meta http-equiv="content-type" content="text/html; charset=utf-8" />
	<link rel="stylesheet" type="text/css" media="screen, projection" href="default.css" />
	<meta name="author" content="Chris Rycroft" />
	<meta name="robots" content="none" />
</head>
<body>
<h1>STZ crack simulation table</h1>
<table>
	<tr id="header">
EOF

@A=split ',',<E>;
$m==0;
foreach (@A) {
	last if $_ eq '';
	$m++;
	s/nu/&nu;/;
	s/varrho/&varrho;/;
	s/mu m/&mu;m/;
	s/_z/<sub>z<\/sub>/;
	s/Gamma/&gamma;/;
	s/Initial chi/&chi;<sub>init<\/sub>/;
	s/ /&nbsp;/g;
	s/GPs/<acronym title="grid points">GP<\/acronym>s/;
	s/mubar/<span style="text-decoration: overline">&mu;<\/span>/;
	print F "\t\t<th>$_</th>\n";
}
print F "\t\t<th><abbr title=\"Minutes\">Min</abbr>/frame</th>\n";
print F "\t\t<th>Frames</th>\n";
print F "\t\t<th>K<sub>I</sub> at p=-2</th>\n";
print F "\t\t<th>Movies</th>\n";
print F "\t</tr>\n";

while(<E>) {
	@A=split ',';
	last if @A[0] eq '';
	print F "\t<tr>\n";
	foreach $i (0..($m-1)) {
		$_=$A[$i];
		s/1\.00E\+0(\d)/10<sup>\1<\/sup>/;
		s/1\.00E-0(\d)/10<sup>-\1<\/sup>/;
		s/(\d)\.00E\+0(\d)/\1 &times; 10<sup>\2<\/sup>/;
		s/(\d)\.00E-0(\d)/\1 &times; 10<sup>-\2<\/sup>/;
		s/(\d)\.50E-0(\d)/\1.5 &times; 10<sup>-\2<\/sup>/;
		s/ /&nbsp;/g;
		s/-/&minus;/g;
		printf F "\t\t<td>$_</td>\n";
	}
	$host='';$cdp='';
	if($_ eq "Derwent") {
		$host="derwent.dhcp.lbl.gov";
	} elsif($_ eq "Yuba") {
		$host="yuba.dhcp.lbl.gov";
	} elsif($_ eq "Merced") {
		$host="avon.dhcp.lbl.gov";
	} elsif($_ eq "Jordan") {
		$host="jordan.lbl.gov";$cdp="/data/chr/";
	} elsif($_ eq "Odra") {
		$host="odra.lbl.gov";
	} elsif($_ eq "Po") {
		$host="po.lbl.gov";
	} elsif($_ eq "Macmini") {
		$host="76.102.87.119";
	}
	if($host ne '') {
		`ssh $host "cd ${cdp}sb/sbsim/$A[0]; ls -l u.*" >temp_timing`;
		open A,"temp_timing";
		$maxf=0;
		$e[$_]=0 foreach 0..100;
		while(<A>) {
			if(m/ (...)  ?(\d*) (\d\d):(\d\d) u\.(\d*)$/) {
				$f=$5;chomp $f;$maxf=$f if $f>$maxf;
				$mo=$ed[$mon{$1}];
				$da=$2;chomp $da;
				$hr=$3;chomp $hr;
				$mi=$4;chomp $mi;
				$e[$f]=$mi+60*($hr+24*($da+$mo));
			}
		}
		close A;
		if($maxf==0) {
			print F "\t\t<td>&#8211;</td>\n\t\t<td>0</td>\n";
		} else {
			$avt=($e[$maxf]-$e[0])/$maxf;
			open A,">temp_timing";
			foreach(1..$maxf) {
				printf A "%g %d\n",$_-0.5,($e[$_]-$e[$_-1]);
			}
			close A;
			open B,">temp_gnuplot";
			print B <<EOF;
set term pngcairo size 800,600
set xlabel "Frame"
set ylabel "Compute time (min)"
unset key
set title "Computation time for $A[0]"
set output '$A[0].png'
set pointsize 1.2
plot [0:*] [0:*] 'temp_timing' w lp pt 7 lt 3
set output
EOF
			close B;
			`gnuplot temp_gnuplot`;
			`mv $A[0].png timing`;
			`chmod a+rx timing/$A[0].png`;
			printf F "\t\t<td><a href=\"timing/$A[0].png\">%.1f</a></td>\n\t\t<td>$maxf</td>\n",$avt;
		}
	} else {
		print F "\t\t<td></td>\n";
	}

	if($host ne '') {
		`rsync -vz $host:${cdp}sb/sbsim/$A[0].fld temp.fld`;
		open B,">temp_gnuplot";
		print B <<EOF;
set term pngcairo size 800,600
set xlabel "K_I (MPa m^0.5)"
set ylabel "Rescaled pressure"
unset key
set title "Minimum pressure for $A[0]"
set output '$A[0]_p.png'
plot [0:80] [-4:0] 'temp.fld' u 2:5 w l lt 3, -2 lt 4
set output
EOF
		close B;
		`gnuplot temp_gnuplot`;
		`mv $A[0]_p.png min_p`;
		`chmod a+rx min_p/$A[0]_p.png`;

		open A,"temp.fld";
		$KI=-1;
		while(<A>) {
			@B=split;
			if($B[4]<-2) {
				$KI=$KK+($B[1]-$KK)*(-2-$pp)/($B[4]-$pp);
				last;
			}
			$pp=$B[4];
			$KK=$B[1];
		}
		close A;
		print F "\t\t<td><a href=\"min_p/$A[0]_p.png\">";
		if($KI<0) {print F "*";}
		else {printf F "%.2f",$KI;}
		print F "</a></td>\n";
	}

	print F "\t\t<td>";
	$first=1;
	foreach $ext ("chi","p","dev") {
		if (-f "$A[0]_$ext.mov") {
			print F "&nbsp;" unless $first;$first=0;
			$ext2=$ext;
			$ext2=~s/chi/&chi;/;
			$ext2=~s/dev/|&sigma;<sub>0<\/sub>|/;
			print F "<a href=\"$A[0]_$ext.mov\">$ext2</a>";
			system "chmod a+rx $A[0]_$ext.mov";
		}
	}
	print F "\t\t</td>\n\t</tr>\n";
}

print F <<EOF;
	</tr>
</table>
</body>
</html>
EOF
