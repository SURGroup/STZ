set pm3d
set size ratio -1
unset surface
set view map
a=4.000001
b=1.000001
xmin=-a;xmax=a;ymin=-b;ymax=b;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x/L'
set ylabel 'y/L'
set term pngcairo size 1100,800
set output 'OUTFILE'
set multiplot
set label 2 at screen 0.148,0.7505 "Direct" tc rgbcolor "#1100ee"
set label 1 at screen 0.148,0.405 "Quasistatic" tc rgbcolor "#1100ee"
set label 3 at screen 0.148,0.075 "Frame FRAME" tc rgbcolor "#880088"
set origin 0,0.35
set size 1,0.5
set cbrange CBRANGE
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
set origin 0,0
set cbrange CBRANGE2
splot [xmin:xmax] [ymin:ymax] 'INFILE2' matrix binary
UND:unset pm3d
UND:set style data lines
UND:unset clabel
UND:set cntrparam levels incremental -4.2,0.4,4.2
UND:set style line 1 linecolor rgb "#CCDDAA" lw 1
UND:set contour
UND:unset surface
UND:set xlabel " "
UND:set ylabel " "
UND:set xtics ""
UND:set ytics ""
UND:unset label 1
UND:unset label 2
UND:unset label 3
UND:set origin 0,0
UND:splot [xmin:xmax] [ymin:ymax] 'XFIELD' matrix binary ls 1
UND:splot [xmin:xmax] [ymin:ymax] 'YFIELD' matrix binary ls 1
UND:set origin 0,0.35
UND:splot [xmin:xmax] [ymin:ymax] 'XFIELD2' matrix binary ls 1
UND:splot [xmin:xmax] [ymin:ymax] 'YFIELD2' matrix binary ls 1
unset multiplot
set output
