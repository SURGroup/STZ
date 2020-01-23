set pm3d
set size ratio -1
unset surface
set view map
a=0.500001
b=1.000001
xmin=-a;xmax=a;ymin=-b;ymax=b;
#xmin=-a*0.5;xmax=a*0.5;
#ymin=-b;
#ymax=0;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 700,900
set output 'OUTFILE'
set multiplot
set label 3 at screen 0.198,0.12 "Frame FRAME" tc rgbcolor "#880088"
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
unset label 3
set style data lines
UND:unset pm3d
UND:unset clabel
UND:set cntrparam levels incremental -2*a,0.2,2*a
UND:set style line 1 linecolor rgb "#CCDDAA" lw 1
UND:set contour
UND:set xlabel " "
UND:set ylabel " "
UND:set xtics ""
UND:set ytics ""
UND:unset label 3
UND:unset surface
UND:splot [xmin:xmax] [ymin:ymax] 'XFIELD' matrix binary ls 1
UND:splot [xmin:xmax] [ymin:ymax] 'YFIELD' matrix binary ls 1
unset surface
set style line 1 linecolor rgb "#FFFFFF" lw 2
set contour
unset clabel
set cntrparam levels discrete 0
unset pm3d
splot [xmin:xmax] [ymin:ymax] 'LEVELSET' matrix binary ls 1
unset multiplot
set output
