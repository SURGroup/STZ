#!/usr/bin/perl
use Getopt::Std;
getopts("dfg:hln:p:tu");

# Print help information if requested
if ($opt_h) {
	print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-f       (Add vector field)\n";
	print "-h       (Print this information)\n";
	print "-l       (Clean the grid up before rendering)\n";
	print "-p <n>   (Render output in parallel using n procs)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-t       (Add tracers)\n";
	print "-u       (Add undeformed fields)\n";
	exit 0;
}

die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

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

# Open rotor mark file if available
if(-e "${e}/rmark") {
	open Q,"${e}/rmark" or die "Error opening mark file\n";
	$rmark=1;
	$rm=-1;
}

# Loop over the available frames
while(-e "${e}/$ARGV[1].$a") {

	# Terminate if the specified frame limit has been reached
	last if defined $opt_n && $a>$opt_n;

	# Prepare output filename
	$za=sprintf "_%04d",$a;
	$of="$e.frames\/fr$za";

	# Get current rotor mark position
	if($rmark) {
		while($rm<$a) {
			($rm,$rmarkx,$rmarky)=split ' ',<Q>;
		}
		die "Mark mismatch\n" unless $rm==$a;
	}

	# Skip existing file if -d option is in use
	if ($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
		print "$a (skipped)\n";
		$a+=$da;
		next;
	}
	print "$a\n";

	# Clean up field if requested
	if($opt_l) {
		system "./clean_output $e $ARGV[1] $a";
		system "mv cleaned_field temp_field$P";
		$infile="temp_field$P";
	} else {
		$infile="$e\/$ARGV[1].$a";
	}

	# Call the utility to clean up the undeformed fields if they are in use
	if($opt_u) {
		system "./clean_output $e X $a -n";
		system "mv cleaned_field temp_X$P";
		$xfield="temp_X$P";
		system "./clean_output $e Y $a -n";
		system "mv cleaned_field temp_Y$P";
		$yfield="temp_Y$P";
	}

	# Prepare vector field and tracer filenames
	$vfield="$e\/vfld.$a";
	$tfile="$e\/trace.$a";

	# Create temporary Gnuplot file
	open C,"template.gnuplot";
	open B,">temp$P.gnuplot";
	while(<C>) {

		# File substitutions for tracers
		if($opt_t) {
			s/TRACERS/$tfile/g;
			s/^TRA://g;
		} else {
			next if m/^TRA:/;
		}

		# File substitutions for undeformed fields
		if($opt_u) {
			s/^UND://;
			s/XFIELD/$xfield/g;
			s/YFIELD/$yfield/g;
		} else {
			next if m/^UND/;
		}

		# File substitutions for fluid vector field
		if($opt_f) {
			s/^VFL://;
			s/VFLD/$vfield/g;
		} else {
			next if m/^VFL/;
		}

		# File substitutions for timing infromation
		if($header) {
			$ti=$t_s+($t_e-$t_s)*$a/$frn;
			$tif=sprintf "%.1f",$ti;
			s/^HEA://;
			s/TIME/$tif/;
		} else {
			next if m/^HEA/;
		}

		# File substitutions for rotor mark
		if($rmark) {
			s/^REC://;
			s/^CIR://;
			s/RMARKX/$rmarkx/;
			s/RMARKY/$rmarky/;
		} else {
			next if m/^REC/;
			next if m/^CIR/;
		}

		# General file substitutions
		s/CBRANGE/$cb/g;
		s/INFILE/$infile/g;
		s/LEVELSET/$e\/phi.$a/g;
		s/OUTFILE/$of/g;
		print B;
	}
	close C;
	close B;

	# Make a fork to create the graph
	exec "gnuplot temp$P.gnuplot" if (($pid[$P]=fork)==0);

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

# Additional code to automatically make a movie
#unless ($opt_m) {
#	print "make_qt 24 $e.frames/fr_0001.png ${e}_$ARGV[1].mov";
#	system "qt_export --sequencerate=20 $e.frames/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${e}_$ARGV[1].mov";
#	system "qt_export --sequencerate=20 $e.frames/fr_0001.png --loadsettings=../../misc/qtprefs/avi --replacefile ${e}_$ARGV[1].avi";
#	system "rm -rf $e.frames";
#}
