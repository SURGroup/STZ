#set palette defined ( 0. "#ff1428", 1/6. "#ff2f92", 1/3. "#ff8ad8", 1/2. "#3791e6", 2/3. "#96e650", 5/6. "#fffb00", 1. "#ffffff")
set size ratio -1
set view map
a=2.00000
b=1.00000
#xmin=-a+3.5;xmax=a+3.5;ymin=-b;ymax=b;
xmin=-a;xmax=a;ymin=-b;ymax=b;
set autoscale zfix
set autoscale cbfix
set view map
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 900,500 font 'Helvetica, 14'
set lmargin at screen 0.1
set rmargin at screen 0.86
set tmargin at screen 0.96
set bmargin at screen 0.12
set output 'OUTFILE'
UND:set multiplot
set label 3 at screen 0.065,0.07 "Frame FRAME" tc rgbcolor "#880088"
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary with image
unset label 3
UND:set style data lines
UND:unset clabel
UND:set style line 1 linecolor rgb "#333333" lw 1
UND:set xlabel " "
UND:set ylabel " "
UND:unset xtics
UND:unset ytics
UND:splot [xmin:xmax] [ymin:ymax] [0:2] 'XFIELD' u 1:2:(1) w l ls 1 
UND:splot [xmin:xmax] [ymin:ymax] [0:2] 'YFIELD' u 1:2:(1) w l ls 1
UND:unset multiplot
set output
