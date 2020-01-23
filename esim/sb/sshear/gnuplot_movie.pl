#!/usr/bin/perl
use Getopt::Std;
getopts("dg:hn:p:");

# Print help information if requested
if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-h       (Print this information)\n";
	print "-p <n>   (Render output in parallel using n procs)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	exit 0;
}

# Check for the correct number of command-line arguments
die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

# Process some options and create output directory
$p=$opt_p?$opt_p:1;$P=0;
$e=@ARGV[0];
mkdir "$e.frames" unless -e "$e.frames";

# If no range is specified, use a default value
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

$a=0;
while(-e "${e}/$ARGV[1].$a") {

	# If the -n option is specified, finish if we have reached the
	# specified number of frames
	if($opt_n) {
		last if $a>$opt_n;
	}

	# Assemble the output filename
	$za=sprintf "_%04d",$a;
	$of="$e.frames\/fr$za";

	# If the -d option is specified, skip files that already exist
	if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
		print "$a (skipped)\n";
		$a++;
		next;
	}

	print "Frame $a\n";

	# Assemble Gnuplot script for this frame
	open C,"template.gnuplot";
	open B,">temp$P.gnuplot";
	$infile="$e\/$ARGV[1].$a";
	while(<C>) {
		s/CBRANGE/$cb/g;
		s/INFILE/$infile/g;
		s/OUTFILE/$of/g;
		print B;
	}
	close C;
	close B;

	# Fork system calls to Gnuplot
	if(fork==0) {
		`gnuplot temp$P.gnuplot`;
		`convert -crop 740x630+100+100 $of.png png24:$of.png`;
		exit;
	}

	# Wait for a batch of forked jobs to finish
	if(++$P==$p) {
		wait foreach 1..$P;$P=0;
	}
	$a++;
}

# Wait for forked jobs to finish
wait foreach 1..$P;

# Commands to make a movie - will not work without additional setup
unless ($opt_f) {
	print "make_qt 15 $e.frames/fr_0001.png ${e}_$ARGV[1].mov";
	system "qt_export --sequencerate=15 $e.frames/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${e}_$ARGV[1].mov";
#	system "rm -rf $e.frames";
}
