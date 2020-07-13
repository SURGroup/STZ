#!/usr/bin/perl
foreach $t (600,630,660) {
	open A,"am${t}l.fin";
	open B,">am${t}n.fin";
	while(<A>) {

#		next if /truncate_model/;
		s/1/0.4/ if /c0/;
#		s/3500/12000/ if /Delta/;
#		s/160/450/ if /Omega/;
#		s/513/1025/ if /gridpoints/;
#		s/9/20/ if /sim_size/;
#		s/150e-6/65e-6/ if /notch_radius/;
#		s/100/160/ if /max_ki/ or /frames/;
#		s/1e-8/5.26038576682752e-9/ if /gamma/;
#		s/2e7/4e7/ if /tmult/;
		print B;
	}
}
