#!/usr/bin/perl

# Load the machines and the corresponding locations of the PSI runs, and save
# the information into a hash table
open A,"places";
while(<A>) {
	($machine,$dir)=split;
	die "Can't process machine address\n" unless $machine=~/^(\w*)\./;
	$ma{$1}=$machine;
	$di{$1}=$dir;
}
close A;

# Load the runs, and sync the key files to the hosts
open B,"jobs";
while(<B>) {
	next if /^( \t)*#/;
	($jn,$ho)=split;
	print "Syncing job $jn to $ho\n";
	die "Hostname $ho not known\n" unless defined $ma{$ho};
	$path="$di{$ho}/psi/run/$jn";
	system "ssh $ma{$ho} mkdir -p $path";
	system "rsync -rz $jn/input?.lmp $jn/master.pov $jn/info.pl $ma{$ho}:$path";
}
