# Crack tip gnuplot file with square view
set pm3d
set size ratio -1
unset surface
set view map
a=5
xmin=-a;xmax=a;ymin=-a;ymax=a;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 720,620
set output 'OUTFILE'
set lmargin at screen 0.05
set rmargin at screen 0.92
set tmargin at screen 0.97
set bmargin at screen 0.12
set multiplot
set label 1 'K_I = FRAME MPa m^{0.5}' left at screen 0.117,0.045 tc lt 1
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
set contour
set cntrparam levels discrete 0
unset pm3d
set style data lines
unset clabel
set style line 1 lw 1 linecolor rgb "#FFFFFF"
set style line 2 linecolor rgb "#CCDDAA"
set pointsize 0.15
splot [xmin:xmax] [ymin:ymax] 'LEVELSET' matrix binary ls 1
UND:set cntrparam levels incremental -a,0.25,a
UND:set style line 1 linecolor rgb "#CCDDAA" lw 1
UND:splot [xmin:xmax] [ymin:ymax] 'XFIELD' matrix binary ls 1
UND:splot [xmin:xmax] [ymin:ymax] 'YFIELD' matrix binary ls 1
unset contour
set surface
unset multiplot
set output
