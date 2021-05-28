CFLAGS += $(HDF5_CFLAGS)
CFLAGS += $(GEOTIFF_CFLAGS)
# Makefile for		: remap
# Module Author		: Orion Lawlor

include ../../make_support/system_rules

OBJS = CLA.o \
	Matrix2D.o \
	mapping.o \
	fetcher.o \
	poly.o \
	quadratic.o \
	mapSamp.o \
	warpMap.o \
	writer.o \
	remap.o

CFLAGS += -Wall $(W_ERROR) $(GEOTIFF_CFLAGS) $(HDF5_CFLAGS)

all: static_library clean

static_library: $(OBJS)
	rm -f libasf_remap.a
	ar r libasf_remap.a $(OBJS)
	$(RANLIB) libasf_remap.a
	mv libasf_remap.a $(LIBDIR)
	cp -p asf_remap.h $(ASF_INCLUDE_DIR)

shared_library: $(OBJS)
	gcc -shared -o libasf_remap.so $(OBJS)
	mv libasf_remap.so $(LIBDIR)
	
clean:
	rm -f *.o
