#! /usr/bin/perl -w
use File::Find;
use File::Path;
$linux=$^O eq "linux";
sub conv;
find {no_chdir=>1,wanted=>\&conv},'.';
print "Making table\n";
#`ssh jordan.lbl.gov "cd public_html/sb; perl make_table.pl; "`;
exit;

sub conv {
	return unless m/\.fld$/;
	s/^.\///;
	s/\.fld//;
	return unless m/^mt/ || m/^iv/ || m/^sc/;
	$fn=$_;
	$opt_u=(m/^sc/)?" -u":"";
	foreach $su ("chi","p","dev") {
		print "$fn $su\n";
		if($linux) {
			`perl crack_movie.pl -n 100 -f -p 3 -c -g 3 -l$opt_u $fn $su`;
			`rsync -rv $fn.frames merced.dhcp.lbl.gov:sb/mov`;
			`ssh merced.dhcp.lbl.gov "cd sb/mov; qt_export --sequencerate=15 $fn.frames/fr_0001.png --loadsettings=../../misc/qtprefs/qt --replacefile ${fn}_$su.mov >/dev/null 2>/dev/null; rm -rf $fn.frames; rsync -rv ${fn}_${su}.mov jordan.lbl.gov:public_html/sb; rm -rf ${fn}_${su}.mov; "`;
		} else {
			`perl crack_movie.pl -n 100 -p 3 -c -g 3 -l$opt_u $fn $su`;
			`rsync -rv ${fn}_${su}.mov jordan.lbl.gov:public_html/sb`;
		}
	}
}
