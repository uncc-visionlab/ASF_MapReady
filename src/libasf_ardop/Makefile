CFLAGS += $(HDF5_CFLAGS)
CFLAGS += $(GEOTIFF_CFLAGS)
include ../../make_support/system_rules

OBJLIB = ardop_libs.o \
	cfft1d_risc.o \
	estdop.o \
	fft_dop.o \
	read_signal.o \
	ardop.o \
	calc_deltas.o \
	get_params.o

SPECLIB = \
	specan.o \
	specan_ml.o \
	specan_patch.o \
	specan_file.o \
	quicklook.o

ASPLIB = patch.o \
	ardop_setup.o \
	rciq.o \
	rmpatch.o \
	acpatch.o \
	calibration.o

# Filter out the ISOC99 FLAGS which causes some problems.
CFLAGS := $(patsubst -D_ISOC99_SOURCE, , $(CFLAGS))
CFLAGS += $(GLIB_CFLAGS) $(W_ERROR) $(JPEG_CFLAGS) $(GEOTIFF_CFLAGS) $(HDF5_CFLAGS)

all: static_library
	mv libasf_ardop.a $(LIBDIR)
	cp -p ardop_defs.h $(ASF_INCLUDE_DIR)
	cp -p read_signal.h $(ASF_INCLUDE_DIR)

static_library: $(OBJLIB) $(ASPLIB)
	make -C ../libasf_export
	rm -f libasf_ardop.a
	ar r libasf_ardop.a $(OBJLIB) $(ASPLIB)
	$(RANLIB) libasf_ardop.a

shared_library: $(OBJLIB) $(ASPLIB)
	$(CC) -shared -o libasf_ardop.so $(OBJLIB) $(ASPLIB)
	mv libasf_ardop.so $(LIBDIR)

clean:
	rm -f core* $(OBJLIB) *.o
