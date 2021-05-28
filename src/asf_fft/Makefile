# Makefile for: fft.a
# Makefile author: Orion Lawlor
# Makefile date: 5/22/98
#

include ../../make_support/system_rules

OBJS =  dxpose.o \
	fft2d.o \
	fftlib.o \
	matlib.o \
	fftext.o 
all: static_library
	
static_library:	$(OBJS)
	rm -f asf_fft.a
	ar rcv asf_fft.a $(OBJS)
	$(RANLIB) asf_fft.a
	mv asf_fft.a $(LIBDIR)
	echo "ASF FFT Library sucessfully built!"
	rm $(OBJS)

shared_library: $(OBJS)
	gcc -shared -o libasf_fft.so $(OBJS)
	mv libasf_fft.so $(LIBDIR)
	
clean:
	-rm -f *.o ../fft.a
