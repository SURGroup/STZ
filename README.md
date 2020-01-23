Eulerian simulation code repository
===================================
This Git repository contains a variety of simulation codes that have been
developed by the Rycroft group @ Harvard and collaborators.

Setup
=====
Most of the projects are written in C++ and can be compiled with the GNU
compiler collection (gcc.gnu.org) using the GNU make utility. The projects
look for a common configuration file called "config.mk" in the top level of
the repository, which sets the compiler name and various compiler options.

To begin, you must create a "config.mk" file from one of the templates
provided in the "config" directory. The simplest method is to create a symlink
(shortcut) from one of the templates. On a Linux system, this can be done by
typing:

ln -s config/config.mk.linux config.mk

On a Mac system, using GCC 9 installed with MacPorts (www.macports.org), use
the following:

ln -s config/config.mk.mac config.mk

For a custom configuration, copy one of the files from the "config" directory
into the root directory, and then manually edit it.

Once the "config.mk" file is created, change into one of the project
directories and type "make" to compile everything.

Related programs
================
No external dependencies are required to compile most of the codes, but
several programs may be useful for analyzing the output.

- Many codes output to formats that can be read by the freeware plotting
  program Gnuplot (www.gnuplot.info). In particular, many output to Gnuplot's
  two-dimensional binary matrix format. Type "help binary matrix" within
  Gnuplot for details on this file format.

- The freeware raytracer POV-Ray (www.povray.org) can be used for high-quality
  renderings for some of the program output files.

- Many of the projects make use of Perl scripts (www.perl.org), for
  post-processing simulation output into movies, and for automatically
  generating complicated source code. Perl has similar functionality to Python.
  It does not have as good numerical libraries as Python, but often gives
  superior performance for text file parsing.

- Some of the programs link to the C++ interface of the PNG library
  (libpng.org) for low-level construction and manipulation of images. libpng
  is easy to install using package management systems on Linux and Mac.

- Many of the C++ files follow the documentation style of the Doxygen sofware
  package (www.doxygen.org), whereby every function within a class is
  documented in a standard comment block beginning with "/**". Doxygen is not
  needed, but this is still a useful documentation standard to follow.

Projects
========

sb - A variety of simulations of plasticity of amorphous materials and bulk
     metallic glasses. Some of these make use of the projection method for
     quasi-static elastoplasticity, originally developed by Rycroft et al.
     "sb" stands for Santa Barbara, where the project started during 2007
     by Rycroft, James Langer, and Frédéric Gibou.

shared - This directory contains common code that is shared among the
         projects. It contains small classes to represent basic vectors and
         matrices in two and three dimensions. It also has a C++ class,
         gp_matrix, for manipulating 2D fields stored in the Gnuplot matrix
         binary format.
         
        
levelset - library for interface tracking using the level set method, using
           the custom algorithms described by Rycroft et al.
           (doi:10.1016/j.jcp.2011.10.009).

tgmg - This is a library to solve large linear systems via a geometric
       multigrid algorithm, making use C++ templates for custom acceleration.
