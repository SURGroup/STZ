#!/usr/bin/perl
use Getopt::Std;
getopts("cdhp:");

if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-c       (Use pngcairo driver instead of postscript)\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-h       (Print this information)\n";
	print "-p <n>   (Render output in parallel using n procs)\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	exit 0;
}

die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

$p=$opt_p?$opt_p:1;$P=0;
$e=@ARGV[0];
mkdir "$e.frames" unless -e "$e.frames";
$cb=$#ARGV==1?"[*:*]":"[@ARGV[2]:@ARGV[3]]";
$a=0;
while(-e "${e}/$ARGV[1].$a") {
	$za=sprintf "_%04d",$a;
	$of="$e.frames\/fr$za";

	if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
		print "$a (skipped)\n";
		$a++;
		next;
	}

	print "$a\n";
	open C,"template.gnuplot";
	open B,">temp$P.gnuplot";
	while(<C>) {
		if ($opt_c) {
			next if m/^EPS:/;
			s/^CAI://g;
		} else {
			next if m/^CAI:/;
			s/^EPS://g;
		}
		s/CBRANGE/$cb/g;
		s/INFILE/$e\/$ARGV[1].$a/g;
		s/LEVELSET/$e\/phi.$a/g;
		s/OUTFILE/$of/g;
		print B;
	}
	close C;
	close B;

	if (fork==0) {
		`gnuplot temp$P.gnuplot`;
		if ($opt_c) {
			`convert -crop 630x540+86+70 $of.png png24:$of.png`;
		} else {
			`nice -n 19 convert -density 500 -crop 1698x1392+464+222 $of.eps png24:$of.png`;
			`nice -n 19 convert -filter lanczos -geometry 630x630 $of.png png24:$of.png`;
			`rm $of.eps`;
		}
		exit;
	}

	if(++$P==$p) {
		wait foreach 1..$P;$P=0;
	}
	$a++;
}
wait foreach 1..$P;
unless ($opt_f) {
	print "make_qt 15 $e.frames/fr_0001.png ${e}_$ARGV[1].mov";
	system "qt_export --sequencerate=15 $e.frames/fr_0001.png --loadsettings=../misc/qtprefs/qt --replacefile ${e}_$ARGV[1].mov";
	system "rm -rf $e.frames";
}
