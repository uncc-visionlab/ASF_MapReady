CFLAGS += $(HDF5_CFLAGS)
CFLAGS += $(GEOTIFF_CFLAGS)
# Makefile for libasf_insar library

include ../../make_support/system_rules

CFLAGS := -Wall $(CFLAGS) $(GSL_CFLAGS) $(PROJ_CFLAGS) $(GLIB_CFLAGS) $(JPEG_CFLAGS) $(GEOTIFF_CFLAGS) $(HDF5_CFLAGS) $(FFTW_CFLAGS) $(XML_CFLAGS)

OBJS = \
	functions.o \
	ips_config.o \
	ips.o \
	complex2polar.o \
	polar2complex.o \
	baseline.o \
	asf_coregister.o \
	asf_igram_coh.o \
	asf_elevation.o \
	asf_phase_unwrap.o \
	asf_baseline.o \
	deramp.o \
	refine_baseline.o \
	phase_filter.o

all: static_library
	mv libasf_insar.a $(LIBDIR)
	cp -p asf_insar.h $(ASF_INCLUDE_DIR)
	cp -p ips.h $(ASF_INCLUDE_DIR)

static_library: $(OBJS)
	rm -f libasf_insar.a
	ar r libasf_insar.a $(OBJS)
	$(RANLIB) libasf_insar.a
	
shared_library: $(OBJS)
	$(CC) -shared -o libasf_insar.so $(OBJS)
	mv libasf_insar.so $(LIBDIR)

$(OBJS): Makefile $(wildcard *.h) $(wildcard ../../include/*.h)

clean:
	rm -rf $(OBJS) libasf_insar.a *~
