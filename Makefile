# Load the common configuration file
include ../config.mk

cextra=$(cflags) --std=c++11

objs=definitions.o
src=$(patsubst %.o,%.cpp,$(objs))
execs=1dfolds 1dfolds_mt

all: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

1dfolds: 1dfolds.cpp $(objs)
	$(cxx) $(cextra) $(gsl_lflags) $(gsl_iflags) -o $@ $^

1dfolds_mt: 1dfolds_mt.cpp $(objs)
	$(cxx) $(cextra) $(gsl_lflags) $(gsl_iflags) -o $@ $^

%.o: %.cpp
	$(cxx) $(cextra) $(gsl_iflags) -c $<

clean:
	rm -f $(objs) $(execs)

.PHONY: clean depend all
