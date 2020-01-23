#!/usr/bin/perl
use Getopt::Std;
getopts("dfg:hlmn:p:tu");

if ($opt_h) {
	print "Usage: ./crack_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
	print "\nOptions:\n";
	print "-d       (Don't duplicate frames that already exist\n";
	print "-f       (Just make the frames, don't make a movie)\n";
	print "-g <n>   (Use an alternative gnuplot header)\n";
	print "-h       (Print this information)\n";
	print "-l       (Don't clean the grid up before rendering)\n";
	print "-n <n>   (Render up to this number of files)\n";
	print "-p <n>   (Override default number of processes to use)\n";
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

# Prepare output directory
$e=@ARGV[0];
mkdir "$e.frames" unless -e "$e.frames";

# Determine range to use
if($#ARGV==3) {
	$cb="[@ARGV[2]:@ARGV[3]]";
} else {
	if($ARGV[1] eq "chi") {
		$cb="[580/21000.:900/21000.]";
		#$cb="[0.025*21/19.0:0.048*21/19.0]";
	} elsif($ARGV[1] eq "dev") {
		$cb="[0:1.4]";
	} elsif($ARGV[1] eq "VS") {
		$cb="[0.995:1.03]"
	} elsif($ARGV[1] eq "p") {
		$cb="[-4.5:0]";
	} else {
		$cb="[*:*]";
	}
}

# Open diagnostic file
open W,"${e}/kfile";
$wa=-1;$a=0;
$P=1;$queue=$nodes==1?1:0;

# Delete old movie if it is present
unlink "${e}_$ARGV[1].mov" unless $opt_m;

while(-e "${e}/$ARGV[1].$a") {

	# Find time and K_I corresponding to the current frame
	do {
		($wa,$wt,$wtt,$wp)=split ' ',<W>;
		last if $wa eq "";
		die "Frame mismatch" if $wa>$a;
	} while($wa<$a);
	last if $wa eq "";
	$tscaled=$wtt<1e-4?($wtt==0?"0 s":sprintf "%.4g µs",$wtt*1e6):
		($wtt<0.1?sprintf "%.4g ms",$wtt*1e3:
			  sprintf "%.4g s",$wtt);
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
	open C,($opt_g?"gp_headers/t$opt_g.gnuplot":"gp_headers/t3.gnuplot");
	open B,">temp$P.gnuplot";
	if($opt_l) {
		$infile="$e\/$ARGV[1].$a";
	} else {
		system "./clean_output $e $ARGV[1] $a";
		system "mv cleaned_field temp_field$P";
		$infile="temp_field$P";
	}
	if($opt_u) {
		system "./clean_output $e X $a";
		system "mv cleaned_field temp_X$P";
		$xfield="temp_X$P";
		system "./clean_output $e Y $a";
		system "mv cleaned_field temp_Y$P";
		$yfield="temp_Y$P";
#		$xfield="$e\/X.$a";
#		$yfield="$e\/Y.$a";
	}
	$lco=$ARGV[1] eq "p"?"#000000":"#ffffff";
	$tfile="$e\/trace.$a";
	while(<C>) {
		next if m/^EPL/;
		if ($opt_t) {
			s/^TRA://g;
			s/TRACERS/$tfile/g;
		} else {
			next if m/^TRA:/;
		}
		next if m/^EPS:/;
		s/^CAI://g;
		if($opt_u) {
			s/^UND://g;
			s/XFIELD/$xfield/g;
			s/YFIELD/$yfield/g;
		} else {
			next if m/^UND/;
		}
		s/CBRANGE/$cb/g;
		s/INFILE/$infile/g;
		s/LEVELSET/$e\/phi.$a/g;
		s/LCOLOR/$lco/;
		s/OUTFILE/$of.png/g;
		s/TIME/$tscaled/g;
		s/SINTEN/$wps/g;
		print B;
	}
	close C;
	close B;

	# Send the POV-Ray file to a node for processing
	exec "gnuplot temp$P.gnuplot; convert -crop 720x620+100+80 $of.png png24:$of.png" if (($pid[$P]=fork)==0);

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
	print "make_qt 15 $e.frames/fr_0001.png ${e}_$ARGV[1].mov";
	system "qt_export --sequencerate=15 $e.frames/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${e}_$ARGV[1].mov";
#	system "rm -rf $e.frames";
}
