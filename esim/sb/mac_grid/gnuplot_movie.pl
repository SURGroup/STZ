#!/usr/bin/perl
use Getopt::Std;
getopts("dg:hln:p:stu");
$gpm="../../shared/gpm_process";

# Print help information if requested
if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-h       (Print this information)\n";
	print "-l       (Clean the grid up before rendering)\n";
	print "-p <n>   (Render output in parallel using n procs)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-s       (Process file into steps)\n";
	print "-t       (Add tracers)\n";
	print "-u       (Add undeformed fields)\n";
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
mkdir "$e.frames" unless -e "$e.frames";

# Set color range
if($#ARGV==3) {
	$cb="[@ARGV[2]:@ARGV[3]]";
} else {
	$cb="[*:*]";
}

$P=1;$a=0;$da=$opt_m?5:1;

# Open header file if available
if(-e "${e}/header") {
	open A,"${e}/header" or die "Error reading header\n";
	$_=<A>;
	($t_s,$t_e,$frn)=split;
	close A;
	$header=1;
}

# Set gnuplot header and read information about image cropping
open A,"gp_headers/t".($opt_g?$opt_g:1).".gnuplot" or die "Can't find gnuplot header";
$_=<A>;
chomp;
$crop=/^# (.*)/?$1:"";

# Read the rest of the file, skipping and altering lines as necessary
$gpn=0;
while(<A>) {

	# Header substitution
	if($header) {s/^HEA://;} else {next if /^HEA:/;}

	# Tracer substitution
	if($opt_t) {s/^TRA://g;} else {next if m/^TRA:/;}
	$gp[$gpn++]=$_;
}
close A;
$gpn--;

# Loop over the available frames
while(-e "${e}/$ARGV[1].$a") {

	# Terminate if the specified frame limit has been reached
	last if defined $opt_n && $a>$opt_n;

	# Prepare output filename
	$za=sprintf "_%04d",$a;
	$of="$e.frames\/fr$za";

	# Skip existing file if -d option is in use
	if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
		print "$a (skipped)\n";
		$a+=$da;
		next;
	}

	# Clean up field if requested
	if(defined $opt_l or defined $opt_s) {
		$infile="temp_field$P";
		system "$gpm ".($opt_l?"-c $e/phi.$a ":"")
		       .($opt_s?"-s ":"")."$e/$ARGV[1].$a $infile";
	} else {
		$infile="$e\/$ARGV[1].$a";
	}

  # Call the utility to clean up the undeformed fields if they are in use
  if($opt_u) {
    system "../../shared/gpm_process -n $e/X.$a temp_X$P";
    $xfield="temp_X$P";
    system "../../shared/gpm_process -n $e/Y.$a temp_Y$P";
    $yfield="temp_Y$P";
  }

	# Prepare tracer filenames
	$tfile="$e\/trace.$a";

	# Create temporary Gnuplot file
	print "Frame $a (thread $P)\n";
	open B,">gp_headers/temp$P.gnuplot";
	foreach $i (0..$gpn) {
		$_=$gp[$i];

		# File substitutions for tracers
		s/TRACERS/$tfile/g if $opt_t;

		# File substitutions for timing information
		if($header) {
			$ti=$t_s+($t_e-$t_s)*$a/$frn;
			$tif=sprintf "%.2f",$ti;
			s/TIME/$tif/;
		}

    # File substitutions for undeformed fields
    if($opt_u) {
      s/^UND://;
      s/XFIELD/$xfield/g;
      s/YFIELD/$yfield/g;
    } else {
      next if m/^UND/;
    }

		# General file substitutions
		s/CBRANGE/$cb/g;
		s/INFILE/$infile/g;
		s/OUTFILE/$of/g;
		print B;
	}
	close B;

	# Make a fork to create the graph
	exec "gnuplot gp_headers/temp$P.gnuplot; convert $crop $of.png png24:$of.png" if (($pid[$P]=fork)==0);

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
wait foreach 1..($queue?$nodes:$h-1);

# Additional code to automatically make a moviedditional code to automatically make a movie
unless ($opt_m) {
  unlink "${e}_$ARGV[1].mov";
  print "make_qt 24 $e.frames/fr_0001.png ${e}_$ARGV[1].mov";
  $uname=`uname`;
  if($uname=~/Linux/) {
    system "ffmpeg -r 24 -i $e.frames/fr_%4d.png -vb 20M ${e}_$ARGV[1].mpg"
  } elsif($uname=~/Darwin/) {
    system "qt_export --sequencerate=24 $e.frames/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${e}_$ARGV[1].mov";
  }
}
