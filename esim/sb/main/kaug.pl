#!/usr/bin/perl
open A,"$ARGV[0]";
open B,">$ARGV[0]_new";

while(<A>) {
	@A=split;
	die "Column count mismatch" unless $#A==1;
	$t=@A[0]*10;
	print B "$A[0] $t $A[1]\n";
}
close A;
close B;
`mv $ARGV[0]_new $ARGV[0]`;
