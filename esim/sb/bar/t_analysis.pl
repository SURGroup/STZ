#!/usr/bin/perl
$|=1;

$ti=0;
while(<>) {
	if(/^(\d*) iters, res (.*)->(.*), (.*) digits per iter$/) {
		($it,$ir,$fr,$digs)=($1,$2,$3,$4);
		$_=<>;
		die "Can't parse input\n" unless /^Sim Time=(.*) Tune=(\d*) Skipped levels=(\d*) Wall clock time=(.*)$/;
		@a=($1,$2,$3,$4);
		$dt=$a[3]-$t;$t=$a[3];
		print "$a[0] $a[1] $a[2] $a[3] $dt $it $ir $fr $digs\n";
	}
}
