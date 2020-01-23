# Include the configuration file
include ../../config.mk

# Lists of files to be built
objs=sbsim.o common.o
src=$(patsubst %.o,%.cc,$(objs))
execs=sshear

all: $(execs)

# Include the file dependencies
include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

sshear: $(objs) sshear.cc
	$(cxx) $(cflags) -o sshear sshear.cc $(objs)

.PHONY: clean depend
