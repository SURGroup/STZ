set pm3d
#set size ratio -1
unset surface
set view map
a=0.500001
b=0.050001
#xmin=-a+3.5;xmax=a+3.5;ymin=-b;ymax=b;
xmin=-a;xmax=a;ymin=-b;ymax=b;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x (m)'
set ylabel 'y (m)'
set cbrange CBRANGE
set term pngcairo size 1200,400
set output 'OUTFILE'
UND:set multiplot
set label 3 at screen 0.148,0.105 "Frame FRAME" tc rgbcolor "#880088"
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
unset label 3
UND:unset pm3d
UND:set style data lines
UND:unset clabel
UND:set cntrparam levels incremental -0.51,0.02,0.51
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
# -crop 990x285+90+70
