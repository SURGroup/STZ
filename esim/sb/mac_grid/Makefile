# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 16th 2014

# Load the common configuration file
include ../../config.mk

objs=common.o mac_sim.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mac_test

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../tgmg lib
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

include Makefile.dep

mac_test: mac_test.cc $(objs)
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
