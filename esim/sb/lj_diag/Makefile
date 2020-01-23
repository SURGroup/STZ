# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : February 10th 2009

# Load the common configuration file
include ../../config.mk

objs=lj_sim.o
src=$(patsubst %.o,%.cc,$(objs))
execs=lj_diag

all: $(execs)

depend:
	$(cxx) -MM $(src) >Makefile.dep

include Makefile.dep

lj_diag: lj_diag.cc lj_sim.o
	$(cxx) $(cflags) -o lj_diag lj_diag.cc lj_sim.o

%.o: %.cc
	$(cxx) $(cflags) -c $<

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all
