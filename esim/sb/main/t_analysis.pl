#!/usr/bin/perl

open A,$ARGV[0] or die "Can't open input file\n";

$n=0;
$ti=0;
$tc=0;
$to=0;

while(<A>) {
	next unless /^# Output frame/;
	@a=split;
	next if $a[3]==0;
	$a[4]=~s/\[//g;
	$a[4]=~s/,//g;
	$ti+=$a[4];
	$tc+=$a[5];
	$to+=$a[7];
	$n+=1;
}

print "Total integration steps: ",$ti;
print "\nTotal computation time: ",$tc;
print "\nTotal time: ",$tc+$to;

print "\n\nAverage computation time/step: ",$tc/$ti;
print "\nAverage output time/step: ",$to/$n,"\n";
