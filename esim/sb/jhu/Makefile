# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : February 10th 2009

# Load the common configuration file
include ../../config.mk

iflags=-I../../tgmg -I../../shared
lflags=-L../../shared -L.

objs=shear_sim.o qs_multi.o common.o bi_interp.o
src=$(patsubst %.o,%.cc,$(objs))
execs=shear_test2 shear_test3 shear_test4 shear_test5 interp_test

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../tgmg lib
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) -I../../shared -I../../tgmg >Makefile.dep

include Makefile.dep

libjhu.a: $(objs)
	rm -f libjhu.a
	ar rs libjhu.a $^

shear_test2: shear_test2.cc libjhu.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -ljhu -lmatvec

shear_test3: shear_test3.cc libjhu.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -ljhu -lmatvec

shear_test4: shear_test4.cc libjhu.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -ljhu -lmatvec

shear_test5: shear_test5.cc libjhu.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -ljhu -lmatvec

interp_test: interp_test.cc
	$(cxx) $(cflags) -o $@ $<

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -I../levelset -c $<

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
