# Load the common configuration file
include ../config.mk

cextra=$(cflags) --std=c++11

objs=definitions.o
src=$(patsubst %.o,%.cpp,$(objs))
execs=1dfolds

all: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

1dfolds: 1dfolds.cpp $(objs)
	$(cxx) $(cextra) -o $@ $^ $(gsl_iflags) $(gsl_lflags)

%.o: %.cpp
	$(cxx) $(cextra) $(gsl_iflags) -c $<

clean:
	rm -f $(objs) $(execs)

.PHONY: clean depend all
