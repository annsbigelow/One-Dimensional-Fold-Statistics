# Load the common configuration file
include ../../config.mk

objs=seg_collect.o
src=$(patsubst %.o,%.cc,$(objs))
execs=sg_test gslrng_test

all: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

sg_test: sg_test.cc $(objs)
	$(cxx) $(cflags) $(gsl_iflags) $(gsl_lflags) -o $@ $^

gslrng_test: gslrng_test.cc
	$(cxx) $(cflags) $(gsl_iflags) $(gsl_lflags) -o $@ $^

%.o: %.cc
	$(cxx) $(cflags) $(gsl_iflags) -c $<

clean:
	rm -f $(objs) $(execs)

.PHONY: clean depend all
