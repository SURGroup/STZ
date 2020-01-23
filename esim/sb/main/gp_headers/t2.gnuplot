#Bar pulling gnuplot file with 3D view
set view 62,189
#xmin=-0.3;xmax=0.3;ymin=0;ymax=0.3;
xmin=-1;xmax=1;ymin=-0.5;ymax=0.5;
#xmin=-0.2;xmax=0.2;ymin=0.05;ymax=0.25;
#xmin=-0.5;xmax=0.5;ymin=-0.25;ymax=0.25;
#xmin=-0.3;xmax=0.3;ymin=-0.15;ymax=0.15;
unset key
set contour
set xlabel 'x'
set ylabel 'y'
set term pngcairo size 1200,850
set output 'OUTFILE.png'
set style data lines
splot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary
set output
