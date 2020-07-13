#!/usr/bin/perl
use Getopt::Std;
getopts("cdfg:hn:p:u");
$gpc="../../shared/make_contour";

# Print help message if -h flag is given
if ($opt_h) {
	print "Usage: ./shear_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "       ./shear_movie.pl -c {<options>} <f1> <s1> <f2> <s2> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	print "-g <n>   (Use an alternative gnuplot header)\n";
	print "-h       (Print this information)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-p <n>   (Override default number of processes to use)\n";
	print "-u       (Add undeformed fields)\n";
	exit 0;
}

# Check command line options
die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3 || defined $opt_c;
die "Need either four or six arguments with -c option" unless $#ARGV==3 || $#ARGV==5 || !defined $opt_c;

# Determine the number of physical cores
if(!defined $opt_p) {
	$uname=`uname`;
	if($uname=~/Linux/) {
		$nodes=`lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l`;
		chomp $nodes;
	} elsif($uname=~/Darwin/) {
		$nodes=`sysctl -n hw.physicalcpu_max`;
		chomp $nodes;
	} else {
		$nodes=4;
	}
} else {
	$nodes=$opt_p;
}

# Prepare output directory
$e=$ARGV[0];
$e2=$ARGV[2] if $opt_c;
$o=($opt_c?"${e}_$e2":$e).".frames";
mkdir $o unless -e $o;

# Determine color bar ranges to use
if($#ARGV==($opt_c?5:3)) {
	$zlo=$opt_c?4:2;
	$zhi=$opt_c?5:3;
	$cb="[$ARGV[$zlo]:$ARGV[$zhi]]";
	$cb2=$cb if $opt_c;
} else {
	$cb=cb_default($ARGV[1]);
	$cb2=cb_default($ARGV[3]) if $opt_c;
}

# Set gnuplot header and link to processing utility
$gph="gp_headers/t".($opt_g?$opt_g:($opt_c?"7":"6")).".gnuplot";
open A,"$gph" or die "Can't find gnuplot header";

# Read the rest of the file, skipping and altering lines as necessary
$gpn=0;
while(<A>) {
	next if m/^EPL/;
	next if m/^EPS:/;
	s/^CAI://g;
	if($opt_u) {
		s/^UND://g;
	} else {
		next if m/^UND/g;
	}
	$gp[$gpn++]=$_;
}
close A;
$gpn--;

# Loop over available files, create a gnuplot script, and render them
$a=0;
$P=0;$queue=$nodes==1?1:0;

while(1) {

	# Check for presence of required files, searching for compressed
	# versions also
	$ex="";
	last if check_file($infile,"$e/$ARGV[1].$a",$ex,$P);
	if($opt_u) {
		last if check_file($xfile,"$e/X.$a",$ex,$P."b");
		last if check_file($yfile,"$e/Y.$a",$ex,$P."c");
	}
	if(defined $opt_c) {
		last if check_file($infile2,"$e2/$ARGV[3].$a",$ex,$P."d");
		if($opt_u) {
			last if check_file($xfile2,"$e2/X.$a",$ex,$P."e");
			last if check_file($yfile2,"$e2/Y.$a",$ex,$P."f");
		}
	}

	# Check for termination condition
	last if defined $opt_n && $a>$opt_n;
	$of=sprintf "$o\/fr_%04d.png",$a;

	# If the -d flag is given, then skip making the image if one already
	# exists
	if ($opt_d && -e $of && -M $infile > -M $of && (!defined $opt_c || -M $infile2 > -M $of)) {
		print "$a (skipped)\n";
		$a++;$infile="$e\/$ARGV[1].$a";
		next;
	}

	# Call the utility to clean up the undeformed fields if they are in use
	if($opt_u) {
		$xfield="$o/temp_X$P";
		$ex.="$gpc $e/X.$a $xfield r -8.2 0.4 42; ";
		$yfield="$o/temp_Y$P";
		$ex.="$gpc $e/Y.$a $yfield r -0.6 0.4 4; ";
		if(defined $opt_c) {
			$xfield2="$o/temp_XX$P";
			$ex.="$gpc $e/X.$a $xfield2 r -8.2 0.4 42; ";
			$yfield2="$o/temp_YY$P";
			$ex.="$gpc $e/Y.$a $yfield2 r -0.6 0.4 4; ";
		}
	}

	# Assemble the gnuplot script for this file
	print "Frame $a (thread $P)\n";
	open B,">$o/temp$P.gnuplot";
	foreach $i (0..$gpn) {
		$_=$gp[$i];
		s/OUTFILE/$of/g;
		s/FRAME/$a/g;
		$stn=sprintf "%.3f",$a*0.025;
		s/STRAIN/$stn/g;

		if($opt_c) {
			if(/CBRANGE2/) {
				print B ($logscale?"":"un")."set logscale cb\n";
				s/CBRANGE2/$cb2/g;
			}
			s/INFILE2/$infile2/g;
			if($opt_u) {
				s/XFIELD2/$xfield2/g;
				s/YFIELD2/$yfield2/g;
			}
		}

		if(/CBRANGE/) {
			print B ($logscale?"":"un")."set logscale cb\n";
			s/CBRANGE/$cb/g;
		}
		s/INFILE/$infile/g;
		if($opt_u) {
			s/XFIELD/$xfield/g;
			s/YFIELD/$yfield/g;
		}
		print B;
	}
	close C;
	close B;

	# Send the POV-Ray file to a node for processing
	exec $ex."gnuplot $o/temp$P.gnuplot 2>/dev/null" if ($pid[$P]=fork)==0;

	# Wait for one of the forked jobs to finish
	if ($queue) {
		$piddone=wait;$P=0;
		$P++ while $piddone!=$pid[$P] && $P<$nodes;
		die "PID return error!" if $P>=$nodes;
	} else {
		$P++;$queue=1 if $P>=$nodes-1;
	}

	# Increment file counter and set new filename
	$a++;
	$infile="$e\/$ARGV[1].$a";
}

# Wait for all remaining forked jobs to complete
wait foreach 1..($queue?$nodes:$P-1);

# Make a movie of the output
unless ($opt_f) {
	$mf="${e}_".($opt_c?($ARGV[1] eq $ARGV[3]?"${e2}_$ARGV[1]":"$ARGV[1]_${e2}_$ARGV[3]"):$ARGV[1]);
	unlink "$mf.mov";
	system "ffmpeg -framerate 30 -y -i $o/fr_%4d.png -preset veryslow -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart $mf.mov";
}

# Delete the temporary output directory
#system "rm -rf $odir";

# Default color bar ranges
sub cb_default {
	if(@_[0] eq "chi") {
		return "[0.028:0.048]";
	} elsif(@_[0] eq "tem") {
		return "[580:850]";
	} elsif(@_[0] eq "dev") {
		return "[0:1.4]";
	} elsif(@_[0] eq "p") {
		return "[-0.3:0.3]";
	} elsif(@_[0] eq "proj") {
		$logscale=1;return "[1e-4:10000]";
	}
	return "[*:*]";
}

# Checks to see if an input file exists. If it doesn't it checks to see if a
# compressed version may exist, in which case it queues up an uncompression
# command.
sub check_file {
	if(-e @_[1]) {@_[0]=@_[1];return 0;}
	@_[0]="temp_mat".@_[3];
	if(-e "@_[1].xz") {
		@_[2].="unxz -c -k @_[1].xz > @_[0]; ";
	} elsif(-e "@_[1].bz2") {
		print "bun\n";
		@_[2].="bunzip2 -c -k @_[1].bz2 > @_[0]; ";
	} elsif(-e "@_[1].gz") {
		print "gun\n";
		@_[2].="gunzip -c @_[1].gz > @_[0]; ";
	} else {
		return 1;
	}
	return 0;
}
