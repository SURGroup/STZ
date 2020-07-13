#!/usr/bin/perl
use Getopt::Std;
getopts("hls");

if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> <frame> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-h       (Print this information)\n";
	print "-l       (Clean the grid up before rendering)\n";
	print "-s       (Stand-alone mode)\n";
	exit 0;
}

die "Need either three or five arguments" unless $#ARGV==2 || $#ARGV==4;

$e=@ARGV[0];
$of="${e}_$ARGV[1]_$ARGV[2]";

if($#ARGV==4) {
	$cb="[@ARGV[3]:@ARGV[4]]";
} else {
	if($ARGV[1] eq "chi") {
		$cb="[0.041:0.045]";
	} elsif($ARGV[1] eq "dev") {
		$cb="[0:1]";
	} else {
		$cb="[-1.8:0]";
	}
}

$mode=$opt_s?"-inc":"";

open C,"bitmap_c.gnuplot" or die "Can't open gnuplot file\n";
open B,">temp.gnuplot";
while(<C>) {
	s/CBRANGE/$cb/;
	s/OUTFILE/${of}_im/;
	print B;
}
close B;
close C;

system "./clean_output $e $ARGV[1] $ARGV[2]";
system "./make_contour $e $ARGV[2]";
system "gnuplot temp.gnuplot";

system "convert -shave 20x20 +matte ${of}_im.png png24:${of}_im.png";

print "pngout\n";
system "pngout ${of}_im.png";
system "convert ${of}_im.png ${of}_im.eps";

print "tempgnu\n";
open C,"gp_headers/t3.gnuplot";
open B,">temp.gnuplot";
while(<C>) {
	next if m/^UND:/;
	next if m/^TRA:/;
	next if m/^EPS:/;
	next if m/^CAI:/;
	next if m/multiplot/;
	next if m/INFILE/;
	s/^EPL://;
	s/ standalone// unless $opt_s;
	s/CBRANGE/$cb/;
	s/LEVELSET/contour/;
	s/OUTFILE/temp/;
	print B;
}
close B;
close C;
system "gnuplot temp.gnuplot";

open C,"temp$mode.eps" or die "Can't open EPS file\n";
while(<C>) {
	print;
	last if m/^\d* \d* N$/;
}
die "Error finding graph edge" unless /^(\d*) (\d*) N$/;
$xlo=$1;
$yhi=$2;
if($opt_s) {
	$_=<C>;die "Error finding graph edge" unless m/^(\d*) (\d*) L$/;$ylo=$2;
} else {
	$_=<C>;die "Error finding graph edge" unless m/^0 -(\d*) V$/;$ylo=$yhi-$1;
}
$_=<C>;die "Error finding graph edge" unless m/^(\d*) 0 V$/;$xhi=$xlo+$1;
close C;
system "mv temp$mode.eps $of$mode.eps";
system "epstopdf $of$mode.eps";

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
	s/{temp$mode}/{$of$mode}/;
	print B;
}
close B;
close C;
