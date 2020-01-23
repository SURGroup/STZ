#!/usr/bin/perl
use Getopt::Std;
getopts("dfhlmn:p:");

if ($opt_h) {
	print "Usage: ./prl_movie.pl {<options>} <filename> <filename2> <suffix>\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	print "-h       (Print this information)\n";
	print "-l       (Don't clean the grid up before rendering)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-p <n>   (Override default number of processes to use)\n";
	exit 0;
}

die "Need three arguments" unless $#ARGV==2;

# Determine the number of processors
if(!defined $opt_p) {
	$uname=`uname`;
	if($uname=~/Linux/) {
		$nodes=`awk '/^processor/ {++n} END {print n+1}' /proc/cpuinfo/`;
		chomp $nodes;
	} elsif($uname=~/Darwin/) {
		$nodes=`sysctl -n hw.ncpu`;
		chomp $nodes;
	} else {
		$nodes=4;
	}
} else {
	$nodes=$opt_p;
}

# Prepare output directory
$e=@ARGV[0];
mkdir "$e.frames" unless -e "$e.frames";

# Open diagnostic file
open W,"${e}/kfile";
$wa=-1;$a=0;
$P=1;$queue=$nodes==1?1:0;

# Delete old movie if it is present
unlink "${e}_$ARGV[2].mov" unless $opt_m;

while(-e "${e}/$ARGV[2].$a") {

	# Find time and K_I corresponding to the current frame
	do {
		($wa,$wt,$wtt,$wp)=split ' ',<W>;
		last if $wa eq "";
		die "Frame mismatch" if $wa>$a;
	} while($wa<$a);
	last if $wa eq "";
	$tscaled=$wtt<1e-4?($wtt==0?"0\\textsf{~s}":sprintf "%.4g\\textsf{~$\mu$s}",$wtt*1e6):
		($wtt<0.1?sprintf "%.4g\\textsf{~ms}",$wtt*1e3:
			  sprintf "%.4g\\textsf{~s}",$wtt);
	$wps=sprintf "%.4g",$wp;
	
	if($opt_n) {
		last if $a>$opt_n;
	}
	$za=sprintf "_%04d",$a;
	$of="$e.frames\/fr$za";

	if ($opt_d && -e $of && -M "$e.$a" > -M $of) {
		print "$a (skipped)\n";
		$a++;
		next;
	}

	#if($a>0) {
#		$pa=$a-1;
#		last if -M "$e\/$ARGV[1].$a" < -M "$e\/$ARGV[1].$pa";
#	}

	print "$a\n";
	open C,"template_prl.gnuplot";
	open B,">temp$P.gnuplot";
	if($opt_l) {
		$infile="$e\/$ARGV[2].$a";
		$infile="$ARGV[1]\/$ARGV[2].$a";
	} else {
		system "./clean_output $e $ARGV[2] $a";
		system "mv cleaned_field temp_field$P";
		$infile="temp_field$P";
		system "./clean_output $ARGV[1] $ARGV[2] $a";
		system "mv cleaned_field temp_fieldb$P";
		$infile2="temp_fieldb$P";
	}
	while(<C>) {
		s/INFILE1/$infile/g;
		s/INFILE2/$infile2/g;
		s/LEVELSET1/$e\/phi.$a/g;
		s/LEVELSET2/$ARGV[1]\/phi.$a/g;
		s/OUTFILE/temp_out$P.tex/g;
		s/TIME/$tscaled/g;
		s/SINTEN/$wps/g;
		print B;
	}
	close C;
	close B;

	# Send the POV-Ray file to a node for processing
	exec "gnuplot temp$P.gnuplot; ./add_key.pl temp_out$P.tex 2; epstopdf temp_out${P}-inc.eps; pdflatex -interaction=batchmode temp_out$P.tex; convert -density 600 -background white -flatten +matte -filter lanczos -geometry 800x800 temp_out$P.pdf png24:$of.png" if (($pid[$P]=fork)==0);

	# Wait for one of the forked jobs to finish
	if ($queue) {
		$piddone=wait;$P=1;
		$P++ while $piddone!=$pid[$P] && $P<=$nodes;
		die "PID return error!\n" if $P>$nodes;
	} else {
		$P++;$queue=1 if $P>=$nodes;
	}
	$a++;
}

wait foreach 1..($queue?$nodes:$h-1);

unless ($opt_f) {
	print "make_qt 15 $e.frames/fr_0001.png ${e}_$ARGV[2].mov";
	system "qt_export --sequencerate=15 $e.frames/fr_0001.png --loadsettings=/Users/chr/misc/qtprefs/qt --replacefile ${e}_$ARGV[2].mov";
#	system "rm -rf $e.frames";
}
