#!/usr/bin/perl
@q=(25,50,75,100,125,150,175,200,225,250);
foreach $t (600,620,640,660) {
	foreach $r (0..$#q) {
		$s=chr 97+$r;
		open A,"cm${t}g.fin";
		open B,">rb${t}$s.fin";
		$ki=50e6*$q[$r]/65.0;
		while(<A>) {
			s/65e-6/$q[$r]e-6/ if /notch_radius/;			
			s/50e6/$ki/ if /typical_ki/;
			s/160/40/ if /frames/;
			s/4/6/ if /tmult/;
			print B;
		}
		close A;
		close B;
	}
}
