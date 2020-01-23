set pm3d
set size ratio -1
unset surface
set view map
a=3.000001
b=1.000001
#xmin=-a+3.5;xmax=a+3.5;ymin=-b;ymax=b;
xmin=-a;xmax=a;ymin=-b;ymax=b;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x (mm)'
set ylabel 'y (mm)'
set cbrange CBRANGE
set term pngcairo size 1100,600 fontscale 1.4
set output 'OUTFILE'
UND:set multiplot
set label 3 at screen 0.115,0.18 "Frame FRAME, strain STRAIN%" tc rgbcolor "#880088"
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
