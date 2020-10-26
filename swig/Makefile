## -*- Makefile -*-
##
## User: arwillis
## Time: Jun 26, 2020 12:34:23 PM
## Makefile created by Oracle Developer Studio.
##
## This file is generated automatically.
##
# 	libatlas-base-dev
# libgdal
CFLAGS += $(HDF5_CFLAGS)
CFLAGS += $(GEOTIFF_CFLAGS)
CFLAGS += $(HDF5_CFLAGS)
# Makefile for libasf_proj library
JAVA_HOME = /usr/lib/jvm/default-java

include ../make_support/system_rules

CFLAGS += $(W_ERROR) \
	$(GLIB_CFLAGS) \
	$(TIFF_CFLAGS) \
	$(GEOTIFF_CFLAGS)

# libs for the test program
LIBS =	$(PROJ_LIBS) \
	$(LIBDIR)/asf.a \
	$(GSL_LIBS) \
	$(TIFF_LIBS) \
	$(GEOTIFF_LIBS) \
	-lm

PYTHON_CFLAGS = $(shell python-config --cflags)
PYTHON_LDFLAGS = $(shell python-config --ldflags)

all: swig
	
## Target: swig
swig: .FORCE
	rm -rf libobjs
	$(MAKE) -C ../src/asf build_only	
	$(MAKE) -C ../src/asf_fft build_only	
	$(MAKE) -C ../src/asf_meta build_only	
	#$(MAKE) -C ../src/asf_view 	
	$(MAKE) -C ../src/libasf_ardop build_only	
	$(MAKE) -C ../src/libasf_convert build_only	
	$(MAKE) -C ../src/libasf_export build_only	
	$(MAKE) -C ../src/libasf_geocode build_only	
	$(MAKE) -C ../src/libasf_import build_only	
	$(MAKE) -C ../src/libasf_insar build_only	
	$(MAKE) -C ../src/libasf_metadata build_only	
	$(MAKE) -C ../src/libasf_proj build_only	
	$(MAKE) -C ../src/libasf_raster build_only	
	$(MAKE) -C ../src/libasf_remap build_only	
	$(MAKE) -C ../src/libasf_sar build_only		
	$(MAKE) -C ../src/libasf_terrcorr build_only	
	$(MAKE) -C ../src/libasf_vector build_only	
	$(MAKE) -C ../src/libifm build_only
	$(MAKE) -C ../src/plan lib
	mkdir -p libobjs
	#cp ../src/asf/*.o libobjs
	#cp ../src/asf_fft/*.o libobjs
	#cp ../src/asf_meta/*.o libobjs	
	#cp ../src/asf_view/asf_view.o ../src/asf_view/read_tiff.o ../src/asf_view/read.o libobjs
	#cp ../src/asf_view/*.o libobjs
	#cp ../src/libasf_ardop/*.o libobjs
	#cp ../src/libasf_convert/*.o libobjs
	#cp ../src/libasf_export/*.o libobjs
	#mv ../src/libasf_import/lut.o ../src/libasf_import/asf_import_lut.o
	#cp ../src/libasf_import/*.o libobjs
	#cp ../src/libasf_insar/*.o libobjs
	#cp ../src/libasf_metadata/*.o libobjs
	#cp ../src/libasf_proj/*.o libobjs
	#cp ../src/libasf_raster/bands.o ../src/libasf_raster/uint8_image.o ../src/libasf_raster/float_image.o ../src/libasf_raster/stats.o libobjs
	#mv ../src/libasf_raster/stats.o ../src/libasf_raster/asf_raster_stats.o
	#mv ../src/libasf_raster/bands.o ../src/libasf_raster/asf_raster_bands.o
	#cp ../src/libasf_raster/*.o libobjs
	#cp ../src/libasf_remap/*.o libobjs
	#cp ../src/libasf_sar/*.o libobjs	
	#cp ../src/libasf_terrcorr/*.o libobjs
	#cp ../src/libasf_vector/kml.o ../src/libasf_vector/utils.o libobjs
	#cp ../src/libifm/dvector.o ../src/libifm/dmatrix.o ../src/libifm/fileExist.o libobjs
	#cp ../src/plan/*.o libobjs
	#rm libobjs/check.o
	#rm libobjs/remap.o
	#rm libobjs/functions.o
	#rm libobjs/phase_filter.o
	#rm libobjs/calibration.o
	$(MAKE) -C ../src/asf clean	
	$(MAKE) -C ../src/asf_fft clean	
	$(MAKE) -C ../src/asf_meta clean
	#$(MAKE) -C ../src/asf_view clean
	$(MAKE) -C ../src/libasf_ardop clean	
	$(MAKE) -C ../src/libasf_convert clean	
	$(MAKE) -C ../src/libasf_export clean	
	$(MAKE) -C ../src/libasf_geocode clean	
	$(MAKE) -C ../src/libasf_import clean	
	$(MAKE) -C ../src/libasf_insar clean	
	$(MAKE) -C ../src/libasf_metadata clean	
	$(MAKE) -C ../src/libasf_proj clean	
	$(MAKE) -C ../src/libasf_raster clean	
	$(MAKE) -C ../src/libasf_remap clean	
	$(MAKE) -C ../src/libasf_sar clean		
	$(MAKE) -C ../src/libasf_terrcorr clean	
	$(MAKE) -C ../src/libasf_vector clean	
	$(MAKE) -C ../src/libifm clean	
	$(MAKE) -C ../src/plan clean	
	mkdir -p ./Python
	echo $(PYTHON_CFLAGS)
	#swig -c++ -python -outdir Python asf_mapready.i
	swig -python -outdir Python asf_mapready.i
	touch Python/__init__.py	
	#gcc -fPIC -c asf_mapready_wrap.cxx $(PYTHON_CFLAGS) -I../include -I/usr/include/geotiff -I/usr/include/hdf5/serial -I/usr/include/libxml2
	gcc -fPIC -c asf_mapready_wrap.c $(PYTHON_CFLAGS) -I../include -I/usr/include/geotiff -I/usr/include/hdf5/serial -I/usr/include/libxml2
	#gcc -fPIC -c asf_mapready_wrap.c -I/usr/include/python2.7 -I../include -I/usr/include/geotiff -I/usr/include/hdf5/serial -I/usr/include/libxml2
	echo $(PYTHON_LDFLAGS)
	#gcc -shared -o Python/_asf_mapready.so asf_mapready_wrap.o $(PYTHON_LDFLAGS) -L../lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lasf -lasf_meta -lasf_import -lasf_insar -lasf_proj -lasf_raster -lhdf5 -lmpi -lgsl -lcblas -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp -lgdal -lglib-2.0 -lstdc++ -lm -lgtk-x11-2.0 -lgdk-x11-2.0 -lpangocairo-1.0 -latk-1.0 -lcairo -lgdk_pixbuf-2.0 -lgio-2.0 -lpangoft2-1.0 -lpango-1.0 -lgobject-2.0 -lglib-2.0 -lfontconfig -lfreetype -lproj
	gcc -shared -o Python/_asf_mapready.so asf_mapready_wrap.o $(PYTHON_LDFLAGS) -Wl,--whole-archive ../lib/asf.a ../lib/asf_meta.a ../lib/libasf_import.a ../lib/libasf_metadata.a ../lib/libasf_proj.a ../lib/libasf_raster.a ../lib/libasf_sar.a  ../lib/asf_fft.a -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp -lcblas -lgdal -lglib-2.0	
	#gcc -shared -o Python/_asf_mapready.so asf_mapready_wrap.o $(PYTHON_LDFLAGS) -Wl,--whole-archive ../lib/asf.a ../lib/asf_meta.a ../lib/libasf_export.a ../lib/libasf_geocode.a ../lib/libasf_import.a ../lib/libasf_insar.a ../lib/libasf_metadata.a ../lib/libasf_proj.a ../lib/libasf_raster.a ../lib/libasf_vector.a ../lib/libifm.a -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lcblas -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp -lgdal -lglib-2.0
	#gcc -shared -o Python/_asf_mapready.so asf_mapready_wrap.o -Wl,--whole-archive ../lib/asf.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/asf_meta.a  -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lcblas -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp -lgdal $(PYTHON_LDFLAGS)
	#gcc -shared -o Python/_asf_mapready.so asf_mapready_wrap.o -Wl,--whole-archive ../lib/asf.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/asf_meta.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_import.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_metadata.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_proj.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_sar.a  -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp 
	mkdir -p ./Java/src/swig/asf
	swig -java -outdir Java/src/swig/asf -package swig.asf asf_mapready.i
	gcc -fPIC -c asf_mapready_wrap.c -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I../include -I/usr/include/geotiff -I/usr/include/hdf5/serial -I/usr/include/libxml2
	#gcc -shared -o Java/libasf_mapready.so asf_mapready_wrap.o -Wl,--whole-archive ../lib/asf.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/asf_meta.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_import.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_proj.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_sar.a  -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp 
	gcc -shared -o Java/src/native_libs/libasf_mapready.so asf_mapready_wrap.o -Wl,--whole-archive ../lib/asf.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/asf_meta.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_import.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_metadata.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_proj.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_sar.a  -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp 
	#gcc -shared -o Python/libasf_mapready.so asf_mapready_wrap.o -Wl,--whole-archive ../lib/asf.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/asf_meta.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_import.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_proj.a  -Wl,--no-whole-archive -Wl,--whole-archive ../lib/libasf_sar.a  -Wl,--no-whole-archive -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/usr/lib/openmpi/lib -lhdf5 -lmpi -lgsl -lxml2 -lproj -ljpeg -ltiff -lgeotiff -lshp 
	#javac Java/swig/asf/*.java
	#javac runme.java -classpath ./Java
	
.FORCE:
	
#### Clean target deletes all generated files ####c
clean:
	rm -rf Java/src/swig/asf/*
	rm -rf Python/asf_mapready
	rm asf_mapready_wrap.o libasf_mapready.so || true



