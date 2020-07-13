#!/usr/bin/perl
open E,"jobs" or die "Can't open job description file\n";
open F,">index.html";

# Pressure thresholds
$pt1=-3.5;
$pt2=-4.5;

print F <<EOF;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN"
"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
	<title>STZ crack simulation table</title>
	<base href="http://math.lbl.gov/~chr/crack/" />
	<meta http-equiv="content-type" content="text/html; charset=utf-8" />
	<link rel="stylesheet" type="text/css" media="screen, projection" href="default.css" />
	<meta name="author" content="Chris Rycroft" />
	<meta name="robots" content="none" />
</head>
<body>
<h1>STZ crack simulation table</h1>
<table>
	<tr>
	<th>Filename</th>
	<th>Model</th>
	<th>T</th>
	<th>&chi;<sub>0</sub></th>
	<th>T<sub>Z</sub></th>
	<th>Grid size</th>
	<th>M</th>
	<th>C. rad. (&mu;m)</th>
	<th>E&nbsp;(GPa)</th>
	<th>&nu;</th>
	<th>s<sub>y</sub>&nbsp;(GPa)</th>
	<th><acronym title="Pulling rate">PR</acronym> (MPa &radic;m/s)</th>
	<th>c<sub>0</sub></th>
	<th>&tau;<sub>0</sub>&nbsp;(s)</th>
	<th>&kappa;</th>
	<th>&Delta;</th>
	<th>&Omega;&nbsp;(&Aring;&sup3;)</th>
	<th>&epsilon;<sub>0</sub></th>
	<th>&chi;<sub>&#8734;</sub></th>
	<th>s<sub>lt</sub></th>
	<th>&Delta;t mult.</th>
	<th>K<sub>I</sub>, p=$pt1</th>
	<th>K<sub>I</sub>, p=$pt2</th>
	<th>Movies</th>
	</tr>
EOF

$q=1;
while(<E>) {

	# Call program to output run setup
	($fn)=split;
	$_=`../crack_digest -h $fn`;
	s/e\+0(\d)/e\+\1/g;
	s/e-0(\d)/e-\1/g;
	s/1e\+(\d*)/10<sup>\1<\/sup>/g;
	s/1e-(\d*)/10<sup>-\1<\/sup>/g;
	s/1\.0*e\+(\d*)/10<sup>\1<\/sup>/g;
	s/1\.0*e-(\d*)/10<sup>-\1<\/sup>/g;
	s/(\d)e\+(\d*)/\1__&times;__10<sup>\2<\/sup>/g;
	s/(\d)e-(\d*)/\1__&times;__10<sup>\2<\/sup>/g;
	s/(\d)\.(\d)e\+(\d*)/\1.\2__&times;__10<sup>\3<\/sup>/g;
	s/(\d)\.(\d)e-(\d*)/\1.\2__&times;__10<sup>\3<\/sup>/g;
	s/__/&nbsp;/g;
	s/-/&minus;/g;
	$fn=~/(\d\d\d)(.*)$/;
	$q^=3 if $l1 ne $1 && $l2 ne $2;
	$l1=$1;$l2=$2;
	$par=$q==2?"odd":"even";
	print F "\t<tr class=\"$par\"><td><a href=\"$fn.fin\">$fn</a></td>$_";
	
	# Create pressure plot	
	open B,">temp_gnuplot";
	print B <<EOF;
set term pngcairo size 800,600
set xlabel "K_I (MPa m^0.5)"
set ylabel "Rescaled pressure"
unset key
set title "Minimum pressure for $fn"
set output 'min_p/${fn}_p.png'
plot [0:100] [-7.5:0] '$fn.fld' u 3:4 w l lt 3, $pt1 lt 4, $pt2 lt 1
set output
EOF
	close B;
	`gnuplot temp_gnuplot`;
	
	# Calculate pressure threshold
	open A,"$fn.fld";
	$KI1=-1;$KI2=-1;$pt1s=0;$pt2s=0;
	while(<A>) {
		@B=split;
		if($B[3]<$pt1&&$pt1s==0) {
			$KI1=$KK+($B[2]-$KK)*($pt1-$pp)/($B[3]-$pp);
			$pt1s=1;
		}
		if($B[3]<$pt2&&$pt2s==0) {
			$KI2=$KK+($B[2]-$KK)*($pt2-$pp)/($B[3]-$pp);
			$pt2s=1;
		}
		last if $pt1s==1 && $pt2s==1;
		$pp=$B[3];
		$KK=$B[2];
	}
	close A;
	print F "<td><a href=\"min_p/${fn}_p.png\">";
	if($KI1<0) {print F "*";}
	else {printf F "%.2f",$KI1;}
	print F "</a></td>";
	print F "<td><a href=\"min_p/${fn}_p.png\">";
	if($KI2<0) {print F "*";}
	else {printf F "%.2f",$KI2;}
	print F "</a></td>";

	# Movies	
	print F "<td>";
	$first=1;
	foreach $exts ("tem","p","dev") {
		$ext=$exts eq "tem" && -f "${fn}_chi.mov"?"chi":$exts;
		if (!-f "${fn}_$ext.mov" && -f "$fn/kfile") {
			#print "$fn $ext\n";
			#`../crack_movie.pl $fn $ext`;
		}
		if (-f "${fn}_$ext.mov") {
			print F "&nbsp;" unless $first;$first=0;
			$ext2=$ext;
			$ext2=~s/chi/&chi;/;
			$ext2=~s/tem/&chi;/;
			$ext2=~s/dev/|&sigma;<sub>0<\/sub>|/;
			$ext2=~s/Dtot/D/;
			$ext2=~s/Dpl/D<sup>pl<\/sup>/;
			$ext2=~s/Del/D<sup>el<\/sup>/;
			print F "<a href=\"${fn}_$ext.mov\">$ext2</a>";
		}
	}
	print F "</td></tr>\n";
}

print F <<EOF;
</table>
</body>
</html>
EOF
