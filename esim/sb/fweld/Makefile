# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 16th 2014

# Load the common configuration file
include ../../config.mk

# Compiler flags for including files and library linking
iflags=-I../../levelset -I../../tgmg -I../../shared
lflags=-L../../levelset -L../../shared -L.

objs=common.o qs_multi.o shear_sim.o
src=$(patsubst %.o,%.cc,$(objs))
execs=shear_test shear_runs shear_interp lin_interp dpl_test

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../levelset
	$(MAKE) -C ../../tgmg lib
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) -I../levelset -Imulti >Makefile.dep

include Makefile.dep

libsbs.a: $(objs)
	rm -f libsbs.a
	ar rs libsbs.a $^

shear_test: shear_test.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_interp: shear_interp.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_runs: shear_runs.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

lin_interp: lin_interp.cc
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags)

dpl_test: dpl_test.cc fweld_model.hh
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
