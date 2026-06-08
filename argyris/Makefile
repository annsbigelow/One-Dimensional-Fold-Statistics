# Include the configuration file
include ../config/config.mk.linux

lflags=`gsl-config --libs`
iflags=`gsl-config --cflags`

# Lists of files to be built
local_objs = mesh.o mesh_io.o pre_comps.o
crumple_objs = rk4.o mesh_en.o common.o ext_potential.o p_grid.o deformation.o
objs = $(local_objs) $(crumple_objs)

src = $(patsubst %.o,%.cc,$(local_objs)) \
      $(patsubst %.o,../crumple/%.cc,$(crumple_objs))

execs=sheet_gen sph_gen unpack mesh_conv convergence argyris_test

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) $(voropp_iflags) $(eigen_iflags) -MM $(src) > Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) libcrump.a $(execs) Makefile.dep

libcrump.a: $(objs)
	rm -f $@
	ar rs $@ $^

$(local_objs): %.o: %.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) $(eigen_iflags) -c $< -o $@

$(crumple_objs): %.o: ../crumple/%.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) $(eigen_iflags) -c $< -o $@

unpack: ../crumple/unpack.cc libcrump.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

sph_gen: ../crumple/sph_gen.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

sheet_gen: ../crumple/sheet_gen.cc
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

mesh_conv: ../crumple/mesh_conv.cc libcrump.a
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

convergence: ../crumple/convergence.cc libcrump.a
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

argyris_test: argyris_test.cc libcrump.a
	$(cxx) $(cflags) $(iflags) $(voropp_iflags) -o $@ $^ $(voropp_lflags) $(lflags)

.PHONY: clean depend
