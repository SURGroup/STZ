# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 16th 2014

# Load the common configuration file
include ../../config.mk

iflags=-I../../levelset -I../../tgmg -I../../shared
lflags=-L../../levelset -L../../shared -L.

objs=common.o qs_multi.o shear_sim.o
src=$(patsubst %.o,%.cc,$(objs))
execs=necking_main

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../levelset
	$(MAKE) -C ../../tgmg lib
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) -I../../levelset -Imulti >Makefile.dep

include Makefile.dep

libsbs.a: $(objs)
	rm -f libsbs.a
	ar rs libsbs.a $^

necking_main: necking_main.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
