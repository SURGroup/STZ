set term pdf color solid
set xlabel 'x'
set ylabel 'y'
f(x)=(x>2?0:(x>1?0.5*(x-1)*(2-x)*(x-2):2*x*x*x-3*x*x+1+0.5*x*x*(1-x)))
g(x)=f(abs(x))
set samples 491
A='cubic_pts'
set pointsize 1.2
set yrange [-0.2:1.2]
set output 'cubic1.pdf'
plot A u 1:2 t 'Data' w p pt 7, g(x) t 'Interpolant' lw 3
set output 'cubic2.pdf'
plot A u 1:3 t 'Data' w p pt 7, g(x)+g(x-1)+g(x-2)+g(x-3)+g(x-4)+g(x-5) t 'Interpolant' lw 3
set output 'cubic3.pdf'
plot A u 1:4 t 'Data' w p pt 7, g(x+4)+g(x+2)+g(x)+g(x-2)+g(x-4) t 'Interpolant' lw 3
set output
