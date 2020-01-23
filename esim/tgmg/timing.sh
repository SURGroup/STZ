for a in 11 13 16 20 40 60 80 120 160 200 240 300 360 420 480 600
do
	g++-mp-4.8 -ansi -pedantic -march=core2 -O3 -fopenmp -DOOO=$a -o poisson poisson.cc
	./poisson | grep TT
done
