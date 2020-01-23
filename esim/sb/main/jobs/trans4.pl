#!/usr/bin/perl
@la=(0,0.2,0.8);
@lb=(0.9,1.4);
@lc=(-4,-4.5);
foreach $t (600,620,630,660) {
	foreach $r (0..11) {
		$s=chr 97+$r;
		$ia=$r%3;
		$ib=($r/3)%2;
		$ic=($r/6)%2;
		open A,"cm${t}g.fin";
		open B,">vo${t}$s.fin";
		while(<A>) {
			s/0.0005/0.05/ if /viscosity/;
			s/20/10/ if /sim_size/;
			s/1025/513/ if /gridpoints/;
			printf B "void_nucl       %g %g 1e-8 %g\n\n",$lc[$ic],$lb[$ib],$la[$ia] if /STZ parameters/;
			print B;
		}
		close A;
		close B;
	}
}
