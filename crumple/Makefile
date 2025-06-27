# Include the configuration file
include ../config.mk

lflags=`gsl-config --libs`
iflags=`gsl-config --cflags`

# Lists of files to be built
objs=rk4.o mesh.o mesh_en.o mesh_io.o common.o ext_potential.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=ode_test ode_conv virt sheet_gen sph_gen unpack mesh_conv ssheet_test bsheet_test

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

libcrump.a: $(objs)
	rm -f $@
	ar rs $@ $^

%.o: %.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -c $<

ode_test: ode_test.cc rk4.o osc.hh
	$(cxx) $(cflags) $(iflags) -o $@ $^

ode_conv: ode_conv.cc rk4.o brusselator.hh
	$(cxx) $(cflags) $(iflags) -o $@ $^

virt: virt.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^

unpack: unpack.cc libcrump.a
	$(cxx) $(cflags) $(iflags) -o $@ $^

sph_gen: sph_gen.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

sheet_gen: sheet_gen.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

mesh_conv: mesh_conv.cc libcrump.a
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

ssheet_test: ssheet_test.cc libcrump.a
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

bsheet_test: bsheet_test.cc libcrump.a
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

.PHONY: clean depend
