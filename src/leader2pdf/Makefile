CFLAGS += $(HDF5_CFLAGS)
CFLAGS += $(GEOTIFF_CFLAGS)
CFLAGS += $(HDF5_CFLAGS)
# Makefile for:   leader2pdf
# Module Author:  Rudi Gens

include ../../make_support/system_rules

LIBS  = \
	$(LIBDIR)/asf_meta.a \
	$(LIBDIR)/asf.a \
	$(GSL_LIBS) \
	-lm

CFLAGS += $(GSL_CFLAGS)

OBJS  = prn_facdr.o \
	prn_raddr.o \
	prn_ifiledr.o \
	prn_mpdr.o \
	prn_dssr.o \
	prn_dhr.o \
	prn_ppdr.o \
	prn_atdr.o \
	prn_rsr.o \
	prn_dqsr.o \
	prn_fdr.o \
	leader2pdf.o

all: metadata
	rm *.o
 
metadata: $(OBJS)
	$(CC) $(CFLAGS) -o leader2pdf $(OBJS) $(LIBS) $(LDFLAGS)
	mv  leader2pdf$(BIN_POSTFIX) $(BINDIR)

clean:
	rm -f *.o core
