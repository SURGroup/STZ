# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 16th 2014

# Load the common configuration file
include ../../config.mk

iflags=-I../../levelset -I../../tgmg -I../../shared
lflags=-L../../levelset -L../../shared -L.

objs=common.o qs_multi.o shear_sim.o bi_interp.o extra.o
src=$(patsubst %.o,%.cc,$(objs))
execs=shear_energy shear_data shear_coord_transform shear_compare shear_energy_adam

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

shear_energy_adam: shear_energy_adam.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_energy: shear_energy.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_data: shear_data.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_coord_transform: shear_coord_transform.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_compare: shear_compare.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

bicubic_vel.cc: bicubic_vel.pl
	./bicubic_vel.pl

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
