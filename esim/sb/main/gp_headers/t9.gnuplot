set pm3d
set size ratio -1
unset surface
set view map
a=2.000001
b=1.000001
xmin=-a;xmax=a;ymin=-b;ymax=b;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 1100,600
set lmargin at screen 0.065
set rmargin at screen 0.86
set tmargin at screen 0.93
set bmargin at screen 0.22
set output 'OUTFILE'
set multiplot
set label 3 at screen 0.065,0.07 "Frame FRAME" tc rgbcolor "#880088"
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
unset label 3
set style data lines
UND:unset pm3d
UND:unset clabel
UND:set cntrparam levels incremental -2*a,0.25,2*a
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
