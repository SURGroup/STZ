# -crop 780x680+50+58
set pm3d
set size ratio -1
unset surface
set view map
a=2
xmin=-a;xmax=a;ymin=-a;ymax=a
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
#set term pngcairo size 860,760
set term pngcairo size 860,760 font "Helvetica, 16"
set output 'OUTFILE.png'
set multiplot
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
unset xtics
unset ytics
unset xlabel
unset ylabel
set contour
set cntrparam levels discrete 0
unset pm3d
set style data lines
unset clabel
set style line 1 lw 2 linecolor rgb "#FFFFFF"
floor=-2+4/255.*49.5
set style line 2 linecolor rgb "#CCDDAA"
set pointsize 0.15
UND:set cntrparam levels incremental -2,0.2,2
UND:set style line 1 linecolor rgb "#CCDDAA" lw 1
UND:splot [xmin:xmax] [ymin:ymax] 'XFIELD' matrix binary ls 1
UND:splot [xmin:xmax] [ymin:ymax] 'YFIELD' matrix binary ls 1
unset contour
HEA:set label 1 't = TIME' left at screen 0.184,0.068 tc rgbcolor "#0000ff" front
set surface
TRA:set style line 2 linecolor rgb "#CCDDAA"
TRA:set pointsize 0.2
TRA:splot [xmin:xmax] [ymin:ymax] 'TRACERS' binary format="%2float" u 1:2:(0) w p ls 2 pt 7 
unset multiplot
set output
