#!/usr/bin/perl

die "Need two arguments" unless $#ARGV==1;

open A,"$ARGV[0]";
open B,">$ARGV[0]-temp";

$f=0;
while(<A>) {
	s/{\\small/[l]{\\small/;
	if(/{\\includegraphics{/) {
		$f=1;
		print B;
		if($ARGV[1]==1) {vertical_bar();}
		elsif($ARGV[1]==2) {vertical_bar2();}
		elsif($ARGV[1]==3) {vertical_bar3();}
		next;
	}
	print B;
}
die "Couldn't modify TeX file\n" if $f==0;
`mv $ARGV[0]-temp $ARGV[0]`;

sub vertical_bar {
	$xl=6300;$xw=180;
	$yl=556;$yw=2776;
	$ts=40;
	$le="-4";
	$l1="-3";
	$l2="-2";
	$l3="-1";
	$l4="0";
	$xtx=$xl+300;

	$xu=$xl+$xw;
	$yu=$yl+$yw;
	$ye=$yl+0.6*$yw/4.6;
	$y1=$yl+1.6*$yw/4.6;
	$y2=$yl+2.6*$yw/4.6;
	$y3=$yl+3.6*$yw/4.6;
	$tsp=$xw+$ts;
	$xwp=0.05*$xw;$ywp=0.05*$yw+0.5;

print B <<EOF;
      \\csname LTb\\endcsname
      \\small
      \\put($xl,$yl){\\includegraphics[width=${xwp}pt,height=${ywp}pt]{colorchartv}}
      \\put($xl,$yl){\\line(1,0){$xw}}
      \\put($xl,$yl){\\line(0,1){$yw}}
      \\put($xl,$yu){\\line(1,0){$tsp}}
      \\put($xu,$yl){\\line(0,1){$yw}}
      \\put($xu,$ye){\\line(1,0){$ts}}
      \\put($xu,$y1){\\line(1,0){$ts}}
      \\put($xu,$y2){\\line(1,0){$ts}}
      \\put($xu,$y3){\\line(1,0){$ts}}
      \\put($xtx,$ye){\\makebox(0,0)[l]{$le}}
      \\put($xtx,$y1){\\makebox(0,0)[l]{$l1}}
      \\put($xtx,$y2){\\makebox(0,0)[l]{$l2}}
      \\put($xtx,$y3){\\makebox(0,0)[l]{$l3}}
      \\put($xtx,$yu){\\makebox(0,0)[l]{$l4}}
EOF


}

sub vertical_bar2 {
	$xl=6300;$xw=180;
	$yl=556;$yw=2776;
	$ts=40;
	$l1="600~K";
	$l2="650~K";
	$l3="700~K";
	$l4="750~K";
	$l5="800~K";
	$l6="850~K";
	$xtx=$xl+280;

	$xu=$xl+$xw;
	$yu=$yl+$yw;
	$y1=$yl+10*$yw/280;
	$y2=$yl+60*$yw/280;
	$y3=$yl+110*$yw/280;
	$y4=$yl+160*$yw/280;
	$y5=$yl+210*$yw/280;
	$y6=$yl+260*$yw/280;
	$tsp=$xw+$ts;
	$xwp=0.05*$xw;$ywp=0.05*$yw+0.5;

print B <<EOF;
      \\csname LTb\\endcsname
      \\small
      \\put($xl,$yl){\\includegraphics[width=${xwp}pt,height=${ywp}pt]{colorchartv}}
      \\put($xl,$yl){\\line(1,0){$xw}}
      \\put($xl,$yl){\\line(0,1){$yw}}
      \\put($xl,$yu){\\line(1,0){$xw}}
      \\put($xu,$yl){\\line(0,1){$yw}}
      \\put($xu,$y1){\\line(1,0){$ts}}
      \\put($xu,$y2){\\line(1,0){$ts}}
      \\put($xu,$y3){\\line(1,0){$ts}}
      \\put($xu,$y4){\\line(1,0){$ts}}
      \\put($xu,$y5){\\line(1,0){$ts}}
      \\put($xu,$y6){\\line(1,0){$ts}}
      \\put($xtx,$y1){\\makebox(0,0)[l]{$l1}}
      \\put($xtx,$y2){\\makebox(0,0)[l]{$l2}}
      \\put($xtx,$y3){\\makebox(0,0)[l]{$l3}}
      \\put($xtx,$y4){\\makebox(0,0)[l]{$l4}}
      \\put($xtx,$y5){\\makebox(0,0)[l]{$l5}}
      \\put($xtx,$y6){\\makebox(0,0)[l]{$l6}}
EOF

}

sub vertical_bar3 {
	$xl=4240;$xw=180;
	$yl=603;$yw=3500;
	$ts=40;
	$l0="-3";
	$l1="-2";
	$l2="-1";
	$l3="0";
	$l4="1";
	$xtx=$xl+300;

	$xu=$xl+$xw;
	$yu=$yl+$yw;
	$y1=$yl+0.25*$yw;
	$y2=$yl+0.5*$yw;
	$y3=$yl+0.75*$yw;
	$tsp=$xw+$ts;
	$xwp=0.05*$xw;$ywp=0.05*$yw;

print B <<EOF;
      \\csname LTb\\endcsname
      \\small
      \\put($xl,$yl){\\includegraphics[width=${xwp}pt,height=${ywp}pt]{colorchartv}}
      \\put($xl,$yl){\\line(1,0){$tsp}}
      \\put($xl,$yl){\\line(0,1){$yw}}
      \\put($xl,$yu){\\line(1,0){$tsp}}
      \\put($xu,$yl){\\line(0,1){$yw}}
      \\put($xu,$y1){\\line(1,0){$ts}}
      \\put($xu,$y2){\\line(1,0){$ts}}
      \\put($xu,$y3){\\line(1,0){$ts}}
      \\put($xtx,$yl){\\makebox(0,0)[l]{$l0}}
      \\put($xtx,$y1){\\makebox(0,0)[l]{$l1}}
      \\put($xtx,$y2){\\makebox(0,0)[l]{$l2}}
      \\put($xtx,$y3){\\makebox(0,0)[l]{$l3}}
      \\put($xtx,$yu){\\makebox(0,0)[l]{$l4}}
EOF

}
