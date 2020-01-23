set pm3d
set size ratio -1
unset surface
set view map
a=4.000001
b=1.000001
#xmin=-a+3.5;xmax=a+3.5;ymin=-b;ymax=b;
xmin=-a;xmax=a;ymin=-b;ymax=b;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 990,290
set lmargin at screen 0.065
set rmargin at screen 0.86
set tmargin at screen 0.93
set bmargin at screen 0.22
set output 'OUTFILE'
UND:set multiplot
set label 3 at screen 0.065,0.07 "Frame FRAME" tc rgbcolor "#880088"
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
unset label 3
UND:unset pm3d
UND:set style data lines
UND:unset clabel
UND:set cntrparam levels incremental -4.2,0.4,4.2
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
UND:unset multiplot
set output
