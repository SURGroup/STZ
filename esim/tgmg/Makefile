# Templated geometric multigrid library (TGMG)
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : December 6th 2015

#Load the common configuration file
include ../config.mk

#List of the common source files
objs=helmholtz.o grating.o poisson_fd.o poisson_fem.o vpoiss_fem.o poisson_1d.o conj_grad.o mpcg_fem.o
src=$(patsubst %.o,%.cc,$(objs))
tgmg_src=tgmg_config.hh tgmg.hh tgmg.cc
execs=pfd_test pfd_time pfem_test pfem_sym basic_cg_test mpcg_test
#helm_test gr_test

#Makefile rules
all: $(execs)

lib: $(tgmg_src)

depend: $(src)
	$(cxx) $(cflags) -MM $(src) -I../levelset -Imulti >Makefile.dep

-include Makefile.dep

gr_test: gr_test.cc grating.o
	$(cxx) $(cflags) -o $@ $^

pfd_test: pfd_test.cc poisson_fd.o
	$(cxx) $(cflags) -o $@ $^

pfd_time: pfd_time.cc poisson_fd.o
	$(cxx) $(cflags) -o $@ $^

pfem_test: pfem_test.cc poisson_fem.o
	$(cxx) $(cflags) -o $@ $^

pfem_sym: pfem_sym.cc poisson_fem.o
	$(cxx) $(cflags) -o $@ $^

helm_test: helm_test.cc helmholtz.o
	$(cxx) $(cflags) -o $@ $^

basic_cg_test: basic_cg_test.cc conj_grad.o poisson_1d.o
	$(cxx) $(cflags) -o $@ $^

mpcg_test: mpcg_test.cc conj_grad.o mpcg_fem.o poisson_fem.o
	$(cxx) $(cflags) -o $@ $^

%.o: %.cc
	$(cxx) $(cflags) -c $<

tgmg.cc: tgmg.in.cc tgmg_compile.pl
	perl tgmg_compile.pl

clean:
	rm -f $(execs) $(objs)

.PHONY: all clean lib depend
