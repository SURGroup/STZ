when gnuplot does not work remember to change the permission to make it executable with chmod

when you want to create a movie of an output, simply run the exec in the main directory and if you want the tracers to be seen, just initialize them in the exec.
if you dont want the tracers (points) to be visible, dont initialize them, and then ignore the errors related to trace when you run the command --->  ./gnuplot_movie.pl -t -u sct_q.out tem

Chris has adjusted the gnuplot_movie.pl to work with the simulations we have. (see GitHub)

the shear_energy_v2.cc has 501 frames instead of the usual 101, that is the only difference
the shear_data_v2.cc is meant to be used to create the movie plots for the displacement

the shear_data_Adam.cc contains the Adam's version of the code (model with less parameters). The shear_energy_Adam.cc is to be used for the EGO code.

shear_energy_Adam2.cc uses an upscale factor to increase the number of steps as a function of the grid resolution.

shear_energy_Adam3.cc is a copy of shear_energy_Adam2.cc but has s_y as a varying parameter

shear_energy_Adam_no_MD uses no b,u0 but initialized the eff. temp. randomly
