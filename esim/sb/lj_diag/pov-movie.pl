#!/usr/bin/perl
use Getopt::Std;
use Sys::Hostname;

getopts("de:hmn:p:rq:s:v");

# Print help message if -h option is specified
if (defined($opt_h)) {
	print "Usage: pov-movie.pl <switches> <snapshot-directory> [<header-number>]\n\n";
	print "Switches:\n";
	print "-d                (Don't duplicate existing files)\n";
	print "-e <num>          (Only render every <num> frame)\n";
	print "-h                (Print this information)\n";
	print "-m                (Automatically create movie)\n";
	print "-n <frame>        (Render up to <frame> frames)\n";
	print "-p <threads>      (Use multiple threads)\n";
	print "-r                (Run threads remotely)\n";
	print "-q <quality>      (Quality of rendering, 1=good, 3=extreme)\n";
	print "-s <frame>        (Render a single frame)\n";
	print "-v                (Verbose output)\n";
	exit 0;
}

die "One or two arguments required" unless @ARGV==1 || @ARGV==2;

$dr=$ARGV[0];
$hn=@ARGV==2?$ARGV[1]:1;
$povm="pov_headers/lj$hn.pov";

push @nlist,"merced.dhcp.lbl.gov" foreach 1..2;
push @nlist,"yuba.dhcp.lbl.gov" foreach 1..4;
push @nlist,"tigris.lbl.gov";
push @nlist,"rhone.lbl.gov" foreach 1..4;
push @nlist,"euphrates.lbl.gov" foreach 1..8;
push @nlist,"klamath.lbl.gov" foreach 1..8;
push @nlist,"jordan.lbl.gov" foreach 1..4;
push @nlist,"snake.lbl.gov" foreach 1..4;
push @nlist,"volga.lbl.gov" foreach 1..3;
push @nlist,"po.lbl.gov" foreach 1..2;

@nlist=("localhost","murray.dhcp.lbl.gov","waveney.dhcp.lbl.gov","derwent.dhcp.lbl.gov");

$host=hostname;
#@nlist=("chr-air.local","chr.local","macmini.local");
@nlist=("imac.local","chr-air.local","chr.local","macmini.local");
#@nlist=("imac.local","chr.local");
$nlist[$_]=~s/$host/localhost/ foreach 0..$#nlist;
$nlist[$_]=~s/imac\.local/10.0.1.8/ foreach 0..$#nlist;
$nlist[$_]=~s/macmini\.local/10.0.1.5/ foreach 0..$#nlist;
$nlist[$_]=~s/chr\.local/10.0.1.6/ foreach 0..$#nlist;
$nlist[$_]=~s/chr-air\.local/10.0.1.9/ foreach 0..$#nlist;
$nlist[$_]=~s/mgaglia\.local/10.0.1.7/ foreach 0..$#nlist;

$nodes=defined $opt_p?$opt_p:1;
$nodes=$#nlist+1 if $opt_p==0;
$queue=$nodes==1?1:0;

$pov_width=580;
$pov_height=580+58;
$pov_flags=$opt_q<=1?($opt_q==1?"+R3 +A0.01 -J":"+R2 +A0.3 -J")
		    :($opt_q==2?"+R5 +A0.01 -J":"+R9 +A0.0001 -J");
#if($opt_q==3) {$pov_width=1400;$pov_height=1400;}
$verb=$opt_v?"":">/dev/null 2>/dev/null";
$every=$opt_e?$opt_e:1;

$h=1;$a=defined $opt_s?$opt_s:0;
$opt_n=$opt_s if defined $opt_s;
while(-e "$dr/g.$a") {

	$fn=sprintf "fr_%04d.png",$a;

	last if defined $opt_n && $a>$opt_n;
	$a+=$every, next if defined $opt_d && -e "$dr/$fn";

	# Assemble the POV file
	open A,">$dr/rtemp$h.pov" or die "Can't open temporary POV file\n";
	open B,$povm or die "Can't open master POV file\n";
	while(<B>) {
		if(/#include "sph\.pov"/) {
			open C,"$dr/g.$a";
			while(<C>) {
				($i,$x,$y,$r,$t)=split;
				$t=5 if $i==12 || $i==10 || $i==109 || $i==66;
				$t=5 if $i==98 || $i==79 || $i==67 || $i==75;
				$t=5 if $i==95 || $i==31 || $i==76 || $i==52;
				$t=5 if $i==59 || $i==90 || $i==9 || $i==64;
				$t=5 if $i==99 || $i==81 || $i==100 || $i==57;
				#$t=3 if $i==50 || $i==35 || $i==107 || $i==34;
				#$t=4 if $i==105|| $i==13 || $i==81 || $i==87;
				printf A "sphere{<%.5f,%.5f,0>,%.5f texture{t$t}}\n",$x,$y,$r*0.9;
				printf A "sphere{<%.5f,%.5f,0>,%.5f texture{t$t}}\n",$x+10,$y,$r*0.9 if $x<-4;
				printf A "sphere{<%.5f,%.5f,0>,%.5f texture{t$t}}\n",$x-10,$y,$r*0.9 if $x>4;
			}
			close C;
			next;
			die;
		}
		print A;
	}
	close B;
	close A;

	# Send the POV-Ray file to a node for processing
	if($opt_r) {
		$hst=$nlist[$h-1];
		print "Frame $a to $hst\n";
		$nice=($hst=~m/localhost/ || $hst=~m/10\.0\.1\./ || $hst=~m/\.local/ || $hst=~m/yuba/ ||
			$hst=~m/po/ || $hst=~m/derwent/ || $hst=~m/merced/ || $hst=~m/waveney/  || $hst=~m/murray/)?"nice -n 19":"nice +19";
		`rsync -rz $dr/rtemp$h.pov $hst:cell/render`;
		exec "ssh $hst \"cd cell/render;$nice povray -D +O$fn +W$pov_width +H$pov_height $pov_flags rtemp$h.pov\" $verb ; rsync -rz $hst:cell/render/$fn $dr ; ssh $hst \"rm cell/render/$fn\" " if (($pid[$h]=fork)==0);
	} else {
		print "Frame $a to fork $h\n";
		exec "cd $dr; povray +O$fn +W$pov_width +H$pov_height $pov_flags rtemp$h.pov $verb" if (($pid[$h]=fork)==0);
	}

	# Wait for one of the forked jobs to finish
	if ($queue) {
		print "Waiting...\n";
		$piddone=wait;$h=1;
		$h++ while $piddone!=$pid[$h]&&$h<=$nodes;
		die "PID return error!\n" if $h>$nodes;
	} else {
		$h++;$queue=1 if $h>=$nodes;
	}

	$a+=$every;
}

wait foreach 1..($queue?$nodes:$h-1);

# Automatically create a movie if requested
`qt_export --sequencerate=40 $dr/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${dr}_c$hn.mov` unless defined $opt_m || defined $opt_s;
