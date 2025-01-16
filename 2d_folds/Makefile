include ../config.mk

iflags=$(gsl_iflags) $(voropp_iflags)
lflags=$(gsl_lflags) $(voropp_lflags)

objs=facet.o sim_flatfold.o mat.o
src=$(patsubst %.o,%.cc,$(objs))
execs=flatfold_test flatfold_gen \
      flatfold_scale flatfold_sc_comb

all:
	$(MAKE) executables

executables: $(execs)

-include Makefile.dep

depend:
	$(cxx) $(cflags) $(iflags) -MM $(src) >Makefile.dep

clean:
	rm -rf $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

flatfold_test: flatfold_test.cc sim_flatfold.o facet.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

flatfold_gen: flatfold_gen.cc sim_flatfold.o facet.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

flatfold_scale: flatfold_scale.cc sim_flatfold.o facet.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

flatfold_sc_comb: flatfold_sc_comb.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

.PHONY: clean depend executables
