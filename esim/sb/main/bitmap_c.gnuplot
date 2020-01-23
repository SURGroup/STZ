set pm3d
set view map
set size ratio -1
unset key
unset colorbox
unset border
unset surface
set rmargin at screen 1
set lmargin at screen 0
set tmargin at screen 1
set bmargin at screen 0
unset xtics
unset ytics
set term pngcairo size 440,440
set cbrange CBRANGE
set output 'OUTFILE.png'
splot [-4.4:4.4] [-4.4:4.4] 'cleaned_field' matrix binary
set output
