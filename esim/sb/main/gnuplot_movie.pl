#!/usr/bin/perl
use Getopt::Std;
getopts("cdg:hln:p:");

if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-c       (Use pngcairo driver instead of postscript)\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-g <n>   (Use an alternative gnuplot header)\n";
	print "-h       (Print this information)\n";
	print "-l       (Clean the grid up before rendering)\n";
	print "-p <n>   (Render output in parallel using n procs)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	exit 0;
}

die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

$p=$opt_p?$opt_p:1;$P=0;
$e=@ARGV[0];
if($ARGV[1] eq "q") {
	$suffix="dev";
	$q_field=1;
} else {
	$suffix=$ARGV[1];
	$q_field=0;
}
mkdir "$e.frames" unless -e "$e.frames";
if($#ARGV==3) {
	$cb="[@ARGV[2]:@ARGV[3]]";
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
$wa=-1;$a=0;

while(-e "${e}/$suffix.$a") {
	do {
		($wa,$wp)=split ' ',<W>;
		die "Wall position mismatch" if $wa>$a;
	} while($wa<$a);

	if($opt_n) {
		last if $a>$opt_n;
	}
	$za=sprintf "_%04d",$a;
	$of="$e.frames\/fr$za";

	if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
		print "$a (skipped)\n";
		$a++;
		next;
	}

	print "$a\n";
	open C,($opt_g?"gp_headers/t$opt_g.gnuplot":"gp_headers/t.gnuplot");
	open B,">temp$P.gnuplot";
	if($opt_l) {
		system "./clean_output $e $suffix $a";
		system "mv cleaned_field temp_field$P";
		$infile="temp_field$P";
	} else {
		$infile="$e\/$suffix.$a";
	}
	while(<C>) {
		next if m/^EPL:/;
		if ($opt_c) {
			next if m/^EPS:/;
			s/^CAI://g;
		} else {
			next if m/^CAI:/;
			s/^EPS://g;
		}
		s/CBRANGE/$cb/g;
		s/'INFILE' /'INFILE' u 1:2:(q(\$3)) /g if $q_field==1;
		s/INFILE/$infile/g;
		s/LEVELSET/$e\/phi.$a/g;
		s/OUTFILE/$of/g;
		s/WALLPOS/$wp/g;
		print B;
	}
	close C;
	close B;

	if (fork==0) {
		`gnuplot temp$P.gnuplot`;
		if ($opt_g==0) {
			if ($opt_c) {
				`convert -crop 885x450+70+70 $of.png png24:$of.png`;
			} else {
				`nice -n 19 convert -density 500 -crop 1698x1392+464+222 $of.eps png24:$of.png`;
				`nice -n 19 convert -filter lanczos -geometry 630x630 $of.png png24:$of.png`;
				`rm $of.eps`;
			}
		} elsif ($opt_g==3) {
			if ($opt_c) {
				`convert -crop 720x620+100+80 $of.png png24:$of.png`;
			}
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
	print "make_qt 15 $e.frames/fr_0001.png ${e}_$suffix.mov";
	system "qt_export --sequencerate=15 $e.frames/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${e}_$ARGV[1].mov";
#	system "rm -rf $e.frames";
}
