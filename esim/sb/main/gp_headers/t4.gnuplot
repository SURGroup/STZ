#Basic gnuplot header file
set pm3d
EPL:set size 1,0.5
set size ratio -1
unset surface
set view map
xmin=-0.75;xmax=0.75;ymin=-0.75;ymax=0.75;
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 960,560
set output 'OUTFILE.png'
set multiplot
q(x)=x>1?(x-1)*(1-1/x):0
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
set style line 2 lw 3 linecolor rgb "#AAFFCC"
set object 1 rect from first -1,-0.5 to first -WALLPOS,0.5
set object 1 rect fc rgb "#101060" fillstyle solid 1.0 lw 0
set object 2 rect from first WALLPOS,-0.5 to first 1,0.5
set object 2 rect fc rgb "#101060" fillstyle solid 1.0 lw 0
aw=0.5
set arrow 1 from first -WALLPOS,ymin to first -WALLPOS,ymax nohead linestyle 2
set arrow 2 from first WALLPOS,ymin to first WALLPOS,ymax nohead linestyle 2
set contour
set cntrparam levels discrete 0
unset pm3d
set style data lines
unset clabel
set style line 1 lw 2 linecolor rgb "#FFFFFF"
splot [xmin:xmax] [ymin:ymax] 'LEVELSET' u 1:2:(abs($1)>WALLPOS?NaN:$3) matrix binary ls 1
unset multiplot
set output
