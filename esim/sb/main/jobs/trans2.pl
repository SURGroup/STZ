#!/usr/bin/perl
@q=(1,1.5,2,4,7,10,15,20,40,70,100,130,160,200);
foreach $t (600,660) {
	foreach $r (0..$#q) {
		$s=chr 97+$r;
		open A,"cm${t}g.fin";
		open B,">pb${t}$s.fin";
		$tm=6e8/$q[$r];
		$tms=sprintf "%.0e",$tm;
		$ga=5.26038576682752e-10*$q[$r];
		while(<A>) {
			s/5.26038576682752e-9/$ga/ if /gamma/;
			s/4e7/$tms/ if /tmult/;
			s/160/40/ if /frames/;
			print B;
		}
		close A;
		close B;
	}
}
