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
#set term pngcairo size 600,192
set term pngcairo size 400,128
set cbrange CBRANGE
set output 'OUTFILE.png'
set object 1 rect from first -1,-0.32 to first -WALLPOS,0.32 front
set object 1 rect fc rgb "#101060" fillstyle solid 1.0 lw 0 front
set object 2 rect from first WALLPOS,-0.32 to first 1,0.32 front
set object 2 rect fc rgb "#101060" fillstyle solid 1.0 lw 0 front
splot [-1:1] [-0.32:0.32] 'cleaned_field' matrix binary
set output
