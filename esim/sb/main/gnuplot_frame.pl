#!/usr/bin/perl
use Getopt::Std;
getopts("cdg:hln:p:");

if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> <frame> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-h       (Print this information)\n";
	print "-l       (Clean the grid up before rendering)\n";
	exit 0;
}

die "Need either three or five arguments" unless $#ARGV==2 || $#ARGV==4;

$e=@ARGV[0];
$of="${e}_$ARGV[1]_$ARGV[2]";

if($#ARGV==4) {
	$cb="[@ARGV[3]:@ARGV[4]]";
} else {
	if($ARGV[1] eq "chi") {
		$cb="[0.07:0.15]";
	} elsif($ARGV[1] eq "dev") {
		$cb="[0:1]";
	} else {
		$cb="[*:*]";
	}
}

open W,"${e}/wall";
$wa=-1;
do {
	($wa,$wp)=split ' ',<W>;
	die "Wall position mismatch" if $wa>$ARGV[2];
} while($wa<$ARGV[2]);
close W;

open C,"bitmap.gnuplot";
open B,">temp.gnuplot";
while(<C>) {
	s/CBRANGE/$cb/;
	s/OUTFILE/${of}_im/;
	s/WALLPOS/$wp/g;print B;
}
close B;
close C;

system "./clean_output $e $ARGV[1] $ARGV[2]";
system "./make_contour $e $ARGV[2]";
system "gnuplot temp.gnuplot";
print "pngout\n";
system "pngout ${of}_im.png";
system "convert ${of}_im.png ${of}_im.eps";

print "tempgnu\n";
open C,"gp_headers/t.gnuplot";
open B,">temp.gnuplot";
while(<C>) {
	next if m/^EPS:/;
	next if m/^CAI:/;
	next if m/multiplot/;
	next if m/INFILE/;
	s/^EPL://;
	s/CBRANGE/$cb/;
	s/LEVELSET/contour/;
	s/OUTFILE/temp/;
	s/WALLPOS/$wp/g;
	print B;
}
close B;
close C;
system "gnuplot temp.gnuplot";

print "slook\n";
open C,"temp.eps";
while(<C>) {
	last if m/^\d* \d* N$/;
}
print "elook\n";
die "Error finding graph edge" unless m/^(\d*) (\d*) N$/;
$xlo=$1;
$yhi=$2;
$_=<C>;die "Error finding graph edge" unless m/^0 -(\d*) V$/;$ylo=$yhi-$1;
$_=<C>;die "Error finding graph edge" unless m/^(\d*) 0 V$/;$xhi=$xlo+$1;
close C;
system "mv temp.eps $of.eps";
system "epstopdf $of.eps";

open C,"temp.tex";
open B,">$of.tex";
while(<C>) {
	$unl=$1 if m/\\setlength{\\unitlength}{(.*)bp}/;
	print B;
	last if m/\\gplbacktext$/;
}
$wid=($xhi-$xlo)*$unl;
$hei=($yhi-$ylo)*$unl;
print B "    \\put($xlo,$ylo){\\includegraphics[width=${wid}bp,height=${hei}bp]{${of}_im}}%\n";
while(<C>) {
	s/{temp}/{$of}/;
	print B;
}
close B;
close C;
