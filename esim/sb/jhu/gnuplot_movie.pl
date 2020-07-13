#!/usr/bin/perl
use Getopt::Std;
getopts("dhj:kn:p:uvw");
$gpm="../../shared/gpm_process";
$gpe="../../shared/gpm_metrics -s";
$gpc="../../shared/make_contour";

# Print help information if requested
if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-d        (Don't duplicate frames that already exist)\n";
	print "-h        (Print this information)\n";
	print "-j <n>    (Render every n frames)\n";
	print "-k        (Keep PNGs and temporary files)\n";
	print "-n <n>    (Render up to this number of files)\n";
	print "-p <n>    (Render output in parallel using n procs)\n";
	print "-u        (Add solid displacement field)\n";
	print "-w        (Render PNGs only without making a movie)\n";
	exit 0;
}

die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

# Determine the number of processors
if(!defined $opt_p) {
	$uname=`uname`;
	if($uname=~/Linux/) {
		$nodes=`awk '/^processor/ {++n} END {print n}' /proc/cpuinfo`;
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

# Make output directory
$e=@ARGV[0];
$ebase=$e;
$ebase=~s/\.out$// or die "Filename should end in '.out'\n";
$odir=$ebase.".frames";
mkdir $odir unless -e $odir;

# Set color range
if($#ARGV==3) {
	$cb="[@ARGV[2]:@ARGV[3]]";
} else {
	$cb="[*:*]";
}

$P=1;$queue=1 if $opt_p==1;$a=0;$k=0;$da=$opt_j?$opt_j:1;

# Open header file if available
if(-e "${e}/header") {
	open A,"${e}/header" or die "Error reading header\n";
	$_=<A>;
	($t_s,$t_e,$frn)=split;
	close A;
	$header=1;
}

# Set gnuplot header and read information about image cropping
open A,"gp_headers/default.gnuplot" or die "Can't find gnuplot header";

# Read the rest of the file, skipping and altering lines as necessary
$gpn=0;
while(<A>) {

	# Header substitution
	if($header) {s/^HEA://;} else {next if /^HEA:/;}

	# Reference map substitution
	if($opt_u) {s/^RMP://;} else {next if /^RMP/;}

	$gp[$gpn++]=$_;
}
close A;
$gpn--;

# Find the grid spacings
($xshift,$yshift)=split ' ',`$gpe ${e}/$ARGV[1].$a`;
$xshift*=0.5;$yshift*=0.5;

# Loop over the available frames
while(-e "${e}/$ARGV[1].$a") {

	# Terminate if the specified frame limit has been reached
	last if defined $opt_n && $a>$opt_n;

	# Prepare output filename
	$za=sprintf "_%04d",$a;
	$of="$odir\/fr$za";

	# Skip existing file if -d option is in use
	if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
		print "$a (skipped)\n";
		$a+=$da;
		next;
	}

	# Compute contours of reference maps if requested
	$ex="";
	if($opt_u) {
		$Xfile="$odir/X$P.dat";
		$ex.="$gpc ${e}/X.$a $Xfile r -4 0.2 41; ";
		$Yfile="$odir/Y$P.dat";
		$ex.="$gpc ${e}/Y.$a $Yfile r -4 0.2 41; ";
	}

	# Create temporary Gnuplot file
	print "Frame $a (thread $P)\n";
	open B,">$odir/temp$P.gnuplot";
	foreach $i (0..$gpn) {
		$_=$gp[$i];

		# File substitutions for timing information
		if($header) {
			$ti=$t_s+($t_e-$t_s)*$a/$frn;
			$tif=sprintf "%.2f",$ti;
			s/TIME/$tif/;
		}

		$infile="$e/$ARGV[1].$a";

		# General file substitutions
		s/XSHIFT/$xshift/g;
		s/YSHIFT/$yshift/g;
		s/CBRANGE/$cb/g;
		s/INFILE/$infile/g;
		s/OUTFILE/$of/g;
		s/XFILE/$Xfile/g;
		s/YFILE/$Yfile/g;
		print B;
	}
	close B;

	# Make a fork to create the graph
	exec $ex."gnuplot -d $odir/temp$P.gnuplot" if (($pid[$P]=fork)==0);

	# Wait for one of the forked jobs to finish
	if ($queue) {
		$piddone=wait;$P=1;
		$P++ while $piddone!=$pid[$P] && $P<=$nodes;
		die "PID return error!\n" if $P>$nodes;
	} else {
		$P++;$queue=1 if $P>=$nodes;
	};
	$a+=$da;
}

# Wait for all the remaining forked jobs to finish
wait foreach 1..($queue?$nodes:$P-1);

# Additional code to automatically make a movie
unless ($opt_w) {
	$mf=$ebase."_".$ARGV[1];
	$uname=`uname`;
	if($uname=~/Linux/) {
		unlink "$mf.mpg";
		system "ffmpeg -r 30 -y -i $odir/fr_%4d.png -vb 20M $mf.mpg"
	} elsif($uname=~/Darwin/) {
		unlink "$mf.mov";
		system "qt_export --sequencerate=30 $odir/fr_0001.png --loadsettings=../../../misc/qtprefs/qt --replacefile $mf.mov";
	}
}

# Delete the temporary output directory
system "rm -rf $odir";
