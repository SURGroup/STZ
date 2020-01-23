# Santa Barbara simulation makefile
#
# Author : Chris H. Rycroft (Harvard SEAS / LBL)
# Email  : chr@alum.mit.edu
# Date   : August 16th 2014

# Load the common configuration file
include ../../config.mk

iflags=-I../../levelset -I../../tgmg -I../../shared
lflags=-L../../levelset -L../../shared -L.

objs=common.o qs_multi.o qs_multi_bd.o bd_sim.o bicubic_vel.o extra_force.o \
     extrap.o fileinfo.o sbfrac.o shear_sim.o dmg_sim.o
src=$(patsubst %.o,%.cc,$(objs))
execs=cavity crack_digest cavity bar_pull wallforce dplas_test \
      instab_test shear_test shear_test2 shear_test3 shear_runs \
      shear_compare shear_mcompare x_sequence bd_test crack_tip y_transect \
      dmg_test shear_random shear_small

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

bar_pull: bar_pull.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o bar_pull bar_pull.cc $(lflags) -lsbs -llevel++

cavity: cavity.cc
	$(cxx) $(cflags) -o $@ $<

cavity2: cavity2.cc stz_model.hh
	$(cxx) $(cflags) -o $@ $<

crack_tip: crack_tip.cc libsbs.a stz_model.hh
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -llevel++

dmg_test: dmg_test.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_test: shear_test.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_small: shear_small.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_random: shear_random.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_runs: shear_runs.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_test2: shear_test2.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_test3: shear_test3.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_compare: shear_compare.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

shear_mcompare: shear_mcompare.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -lmatvec

bd_test: bd_test.cc libsbs.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lsbs -llevel++

instab_test: instab_test.cc libsbs.a stz_model.hh
	$(cxx) $(cflags) -o $@ $< $(lflags) -lsbs

crack_digest: crack_digest.cc libsbs.a stz_model.hh
	$(cxx) $(cflags) -o $@ $< $(lflags) -lsbs

dplas_test: dplas_test.cc libsbs.a stz_model.hh
	$(cxx) $(cflags) -o $@ $< $(lflags) -lsbs

y_transect: y_transect.cc
	$(cxx) $(cflags) $(iflags) $(png_iflags) -o $@ $< $(lflags) -lgpmtx $(png_lflags)

x_sequence: x_sequence.cc
	$(cxx) $(cflags) $(iflags) $(png_iflags) -o $@ $< $(lflags) -lgpmtx $(png_lflags)

noise_pull: noise_pull.cc $(SOURCE)
	$(cxx) $(cflags) -o $@ $<

thread_test: thread_test.cc
	$(cxx) $(cflags) -o $@ $< -lpthread

transect: transect.cc
	$(cxx) $(cflags) -o $@ $<

wallforce: wallforce.cc
	$(cxx) $(cflags) -o $@ $<

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

bicubic_vel.cc: bicubic_vel.pl
	./bicubic_vel.pl

clean:
	rm -f $(execs) $(objs)

.PHONY: clean all executables depend
