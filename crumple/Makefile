# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=rk4.o mesh.o mesh_en.o mesh_io.o common.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=ode_test ode_conv virt sheet_gen # sph_gen mesh_conv ssheet_test bsheet_test

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
	$(cxx) $(cflags) $(voropp_iflags) -c $<

ode_test: ode_test.cc rk4.o osc.hh
	$(cxx) $(cflags) -o $@ $^

ode_conv: ode_conv.cc rk4.o brusselator.hh
	$(cxx) $(cflags) -o $@ $^

virt: virt.cc
	$(cxx) $(cflags) -o $@ $^

unpack: unpack.cc libcrump.a
	$(cxx) $(cflags) -o $@ $^

sph_gen: sph_gen.cc
	$(cxx) $(cflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags)

sheet_gen: sheet_gen.cc
	$(cxx) $(cflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags)

mesh_conv: mesh_conv.cc libcrump.a
	$(cxx) $(cflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags)

ssheet_test: ssheet_test.cc libcrump.a
	$(cxx) $(cflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags)

bsheet_test: bsheet_test.cc libcrump.a
	$(cxx) $(cflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags)

.PHONY: clean depend
