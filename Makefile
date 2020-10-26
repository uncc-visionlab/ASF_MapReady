# This make file lists all the system components used by the latest
# tools in one flat list, making no attempt to sort out the build
# order automatically.

# Note that the order in which the modules are listed is critical!
# Note also that this simplistic build scheme takes the traditional
# built-too-much approach: building the default target all will
# rebuild everything.  So if you want to be sure all the prerequisite
# to your tool are rebuilt, you should do a make all.  On the other
# hand, if you want to rebuild just the module you are working on,
# just do a make in that subdirectory.

MAPREADY_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/asf_fft \
	src/libasf_raster \
	src/libasf_sar \
	src/libasf_import \
	src/libasf_vector \
	src/libasf_geocode \
	src/libasf_export \
	src/libasf_ardop \
	src/libasf_terrcorr \
	src/libasf_insar \
	src/libasf_remap \
	src/asf_import \
	src/asf_export \
	src/asf_geocode \
	src/asf_terrcorr \
	src/asf_calpol \
	src/asf_calibrate \
	src/asf_gamma_import \
	src/asf_airsar_import \
	src/asf_polsarpro_import \
	src/libasf_convert \
	src/asf_mapready \
	src/refine_geolocation \
	src/shift_geolocation \
	src/flip \
	src/write_ppf \
	src/analyze_yaw \
	src/adjust_bands \
	src/fftMatch \
	src/trim \
	src/asf_subset \
	src/sr2gr \
	src/gr2sr \
	src/to_sr \
	src/deskew \
	src/libasf_metadata \
	src/metadata \
	src/resample \
	src/fill_holes \
	src/meta2envi \
	src/meta2xml \
	src/mosaic \
	src/llh2ls \
	src/smooth \
	src/farcorr \
	src/diffmeta \
	src/diffimage \
	src/geoid_adjust \
	src/sgpsdp \
	src/plan \
	src/brs2jpg \
	src/add_aux_band \
	src/fit_warp \
	src/remap \
	src/update_state \
	src/clm \
	src/byteswap \
	src/write_hdf5_xml \
	src/stats \
	src/asf2geobrowse \
	src/sqrt_img \
	src/color_browse \
	src/color_diff \
	src/cleanup_pixel \
	src/annotate_image \
	src/populate_meta_field \
	src/trim_wedges \
	src/kernel \
	src/make_overlay \
	src/asf_kml_overlay \
	src/sample_plugin \
	src/convert2vector \
	src/mapready src/metadata_gui src/asf_view src/proj2proj

MAPREADY_SOURCE = \
	$(MAPREADY_MODULES) \
	include \
	doc/mapready_manual.pdf \
	doc/mapready-version_history.txt \
	doc/arcgis_mask.pdf \
	doc/dem_manual.pdf \
	COPYING \
	make_support \
	configure \
	README_src.txt \
	Makefile.in

STP_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/asf_fft \
	src/libasf_raster \
	src/libasf_sar \
	src/libasf_import \
	src/libasf_vector \
	src/libasf_geocode \
	src/libasf_export \
	src/libasf_ardop \
	src/libasf_terrcorr \
	src/libasf_insar \
	src/metadata \
	src/asf_import \
	src/asf_export \
	src/sgpsdp \
	src/plan \
	src/asf_view \
	src/metadata \
	src/metadata_gui \
	src/stp

AP_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/asf_fft \
	src/libasf_raster \
	src/libasf_sar \
	src/libasf_import \
	src/libasf_vector \
	src/libasf_geocode \
	src/libasf_export \
	src/libasf_terrcorr \
	src/libasf_ardop \
	src/libasf_convert \
	src/metadata \
	src/metadata_gui \
	src/sgpsdp \
	src/plan \
	src/asf_view

REQ_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/libasf_raster \
	src/metadata \
	src/req

C2V_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/libasf_raster \
	src/libasf_sar \
	src/libasf_import \
	src/libasf_vector \
	src/libasf_geocode \
	src/libasf_export \
	src/convert2vector \
	src/convert2vector_gui

JPL_MODULES = \
	contrib/jplmosaicsuite \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/byteswap \
	src/chgndv

CREATE_THUMBS_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/asf_fft \
	src/libasf_raster \
	src/libasf_sar \
	src/libasf_import \
	src/libasf_vector \
	src/libasf_export \
	src/libasf_ardop \
	src/create_thumbs

DATA_QC_MODULES = \
	src/asf \
	src/libasf_proj \
	src/asf_meta \
	src/data_qc

FARADAY_MODULES = \
	src/asf \
	src/faraday_prediction

IGNORED_INCLUDES = \
	include/ardop_defs.h \
	include/asf_baseline.h \
	include/asf_complex.h \
	include/asf_convert.h \
	include/asf_export.h \
	include/asf_geocode.h \
	include/asf_glib.h \
	include/asf.h \
	include/asf_import.h \
	include/asf_insar.h \
	include/asf_jpeg.h \
	include/asf_meta.h \
	include/asf_raster.h \
	include/asf_sar.h \
	include/asf_simulation.h \
	include/asf_terrcorr.h \
	include/asf_tiff.h \
	include/asf_vector.h \
	include/banded_float_image.h \
	include/calibrate.h \
	include/caplib.h \
	include/date.h \
	include/dateUtil.h \
	include/doppler.h \
	include/envi.h \
	include/doppler.h \
	include/float_image.h \
	include/frame_calc.h \
	include/gamma.h \
	include/get_ceos_names.h \
	include/get_stf_names.h \
	include/ips.h \
	include/libasf_proj.h \
	include/line_header.h \
	include/matrix.h \
	include/metadisplay.h \
	include/meta_init_stVec.h \
	include/meta_project.h \
	include/plan.h \
	include/plan_internal.h \
	include/polygon.h \
	include/poly.h \
	include/read_signal.h \
	include/sgpsdp.h \
	include/spheroids.h \
	include/uint8_image.h \
	include/ursa.h \
	include/vector.h \
	include/xml_util.h

S_LIBDIR   = lib
S_BINDIR   = bin
S_SHAREDIR = share/asf_tools
S_DOCDIR   = $(S_SHAREDIR)/doc

ALL_MODULES = \
	$(MAPREADY_MODULES) \
	$(STP_MODULES) \
	$(AP_MODULES) \
	$(REQ_MODULES) \
	$(C2V_MODULES) \
	$(CREATE_THUMBS_MODULES) \
	$(DATA_QC_MODULES)

MAPREADY_SHAREDIR = $(S_SHAREDIR)/mapready
STP_SHAREDIR = $(S_SHAREDIR)/stp
AP_SHAREDIR = $(S_SHAREDIR)/ap
REQ_SHAREDIR = $(S_SHAREDIR)/req
C2V_SHAREDIR = $(S_SHAREDIR)/c2v
JPL_SHAREDIR = $(S_SHAREDIR)/jpl
CREATE_THUMBS_SHAREDIR = $(S_SHAREDIR)/create_thumbs
DATA_QC_SHAREDIR = $(S_SHAREDIR)/data_qc

VERSION_MAPREADY = $(shell awk -F '"' '$$1 ~ /define MAPREADY_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_ASPS = $(shell awk -F '"' '$$1 ~ /define ASPS_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_STP = $(shell awk -F '"' '$$1 ~ /define STP_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_AP = $(shell awk -F '"' '$$1 ~ /define AP_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_REQ = $(shell awk -F '"' '$$1 ~ /define REQ_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_C2V = $(shell awk -F '"' '$$1 ~ /define C2V_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_JPL = $(shell awk -F '"' '$$1 ~ /define JPL_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_CREATE_THUMBS = $(shell awk -F '"' '$$1 ~ /define CREATE_THUMBS_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_DATA_QC = $(shell awk -F '"' '$$1 ~ /define DATA_QC_VERSION_STRING/ {print $$2}' include/asf_version.h)
VERSION_FARADAY = $(shell awk -F '"' '$$1 ~ /define FARADAY_STRING/ {print $$2}' include/asf_version.h)

# Build the MapReady package by default
default: mapready

# MapReady Package
mapready: mkdirs_for_build
	mkdir -p $(MAPREADY_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_MAPREADY > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(MAPREADY_SHAREDIR) >> make_support/build_pkg.mk
	cp doc/mapready_manual.pdf $(S_DOCDIR)
	cp doc/mapready-version_history.txt $(S_DOCDIR)
	cp doc/arcgis_mask.pdf $(S_DOCDIR)
	cp doc/dem_manual.pdf $(S_DOCDIR)
	cp COPYING $(MAPREADY_SHAREDIR)
	echo $(VERSION_MAPREADY) > $(MAPREADY_SHAREDIR)/mapready_version.txt
	$(foreach MODULE, $(MAPREADY_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX ASF MapReady Package Compiled! XXXXXXXXXXX"
	@ echo ""

# MapReady Package for SEASAT processor
asps: mkdirs_for_build
	mkdir -p $(MAPREADY_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_ASPS > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(MAPREADY_SHAREDIR) >> make_support/build_pkg.mk
	cp doc/mapready_manual.pdf $(S_DOCDIR)
	cp doc/mapready-version_history.txt $(S_DOCDIR)
	cp doc/arcgis_mask.pdf $(S_DOCDIR)
	cp doc/dem_manual.pdf $(S_DOCDIR)
	cp COPYING $(MAPREADY_SHAREDIR)
	echo $(VERSION_ASPS) > $(MAPREADY_SHAREDIR)/mapready_version.txt
	$(foreach MODULE, $(MAPREADY_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX ASF ASPS Package Compiled! XXXXXXXXXXX"
	@ echo ""

mapready_source:
	mkdir asf_tools
# This is essentially a cp --parents but not all platforms have the --parents option
	tar --exclude=.svn -cf - $(MAPREADY_SOURCE) | tar xf - -C asf_tools
	tar -czvf mapready-src.tar.gz asf_tools
	rm -fr asf_tools

# SAR Training Processor Package
stp: mkdirs_for_build
	mkdir -p $(STP_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_STP > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(STP_SHAREDIR) >> make_support/build_pkg.mk
	cp doc/SAR_training_processor.pdf $(S_DOCDIR)
	cp doc/stp-version_history.txt $(S_DOCDIR)
	cp COPYING $(STP_SHAREDIR)
	echo $(VERSION_STP) > $(STP_SHAREDIR)/stp_version.txt
	$(foreach MODULE, $(STP_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX ASF SAR Training Processor Package Compiled! XXXXXXXXXXX"
	@ echo ""

# Acquisition Request Planner Package
ap: mkdirs_for_build
	mkdir -p $(AP_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_AP > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(AP_SHAREDIR) >> make_support/build_pkg.mk
	cp COPYING $(AP_SHAREDIR)
	echo $(VERSION_AP) > $(AP_SHAREDIR)/acquisitionplanner_version.txt
	$(foreach MODULE, $(AP_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX ASF Acquisition Planning Package Compiled! XXXXXXXXXXX"
	@ echo ""

# Request Generator Package
req: mkdirs_for_build
	mkdir -p $(REQ_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_REQ > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(REQ_SHAREDIR) >> make_support/build_pkg.mk
	cp COPYING $(REQ_SHAREDIR)
	echo $(VERSION_REQ) > $(REQ_SHAREDIR)/req_version.txt
	$(foreach MODULE, $(REQ_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX ASF Request Generator Compiled! XXXXXXXXXXX"
	@ echo ""

# Convert To Vector Package
c2v: mkdirs_for_build
	mkdir -p $(C2V_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_C2V > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(C2V_SHAREDIR) >> make_support/build_pkg.mk
	cp doc/ConvertToVector.pdf $(S_DOCDIR)
	cp doc/ConvertToVector-version_history.txt $(S_DOCDIR)
	cp COPYING $(C2V_SHAREDIR)
	echo $(VERSION_C2V) > $(C2V_SHAREDIR)/c2v_version.txt
	$(foreach MODULE, $(C2V_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX Convert To Vector Compiled! XXXXXXXXXXX"
	@ echo ""

# JPL Mosaic Suite
jpl: mkdirs_for_build
	mkdir -p $(JPL_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_JPL> make_support/build_pkg.mk
	echo SHAREDIR = ../../$(JPL_SHAREDIR) >> make_support/build_pkg.mk
	cp COPYING $(JPL_SHAREDIR)
	echo $(VERSION_JPL) > $(JPL_SHAREDIR)/jplmosaicsuite_version.txt
	$(foreach MODULE, $(JPL_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX JPL Mosaic Suite Compiled! XXXXXXXXXXX"
	@ echo ""

# Create thumbs
create_thumbs: mkdirs_for_build
	mkdir -p $(CREATE_THUMBS_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_CREATE_THUMBS > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(CREATE_THUMBS_SHAREDIR) >> make_support/build_pkg.mk
	cp doc/create_thumbs-version_history.txt $(S_DOCDIR)
	cp COPYING $(CREATE_THUMBS_SHAREDIR)
	echo $(VERSION_CREATE_THUMBS) > $(CREATE_THUMBS_SHAREDIR)/create_thumbs_version.txt
	$(foreach MODULE, $(CREATE_THUMBS_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX Create Thumbs package Compiled! XXXXXXXXXXX"
	@ echo ""

# Data QC
data_qc: mkdirs_for_build
	mkdir -p $(DATA_QC_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_DATA_QC > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(DATA_QC_SHAREDIR) >> make_support/build_pkg.mk
	cp COPYING $(DATA_QC_SHAREDIR)
	echo $(VERSION_DATA_QC) > $(DATA_QC_SHAREDIR)/data_qc_version.txt
	$(foreach MODULE, $(DATA_QC_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX Data QC Compiled! XXXXXXXXXXX"
	@ echo ""

# Faraday prediction
faraday: mkdirs_for_build
	mkdir -p $(FARADAY_SHAREDIR)
	-chmod -R ug+w share
	echo BUILD_PKG = -D_PKG_FARADAY > make_support/build_pkg.mk
	echo SHAREDIR = ../../$(FARADAY_SHAREDIR) >> make_support/build_pkg.mk
	cp doc/faraday_prediction-version_history.txt $(S_DOCDIR)
	cp COPYING $(FARADAY_SHAREDIR)
	echo $(VERSION_FARADAY) > $(FARADAY_SHAREDIR)/faraday_version.txt
	$(foreach MODULE, $(FARADAY_MODULES), $(MAKE) -C $(MODULE) &&) true
	@ echo ""
	@ echo "  XXXXXXXXXXXX Faraday prediction Compiled! XXXXXXXXXXX"
	@ echo ""

# Build using our old system
oldtools: mkdirs_for_build
	cd make_support; $(MAKE); ./makemake linux; cd ..
	$(MAKE) --makefile=Makefile.old

# Get a nice version number for our package.
DEBIAN_VERSION = $(VERSION_MAPREADY)
ifeq ($(DEBIAN_VERSION),MAKE-DEV)
	DEBIAN_VERSION = 0~$$(date +%Y%m%d)
endif
DATEFORMAT = "%a, %d %b %Y %T %z"

# Generate a debian package.
deb: default
	if [ -d asf-mapready ]; then \
		rm -r asf-mapready; \
	fi
	mkdir -p asf-mapready/DEBIAN
	$(MAKE) DESTDIR=asf-mapready install
	mkdir -p "asf-mapready/usr/share/doc/asf-mapready"
	# Generate the changelog. I'm not looking forward to this.
	echo "asf-mapready ($(DEBIAN_VERSION)) unstable; urgency=low" \
		> asf-mapready/usr/share/doc/asf-mapready/changelog
	git log -n 10 | \
		sed -e "s/^/    /" \
		    -e "/^    commit/s/^    /  * /" \
		>> asf-mapready/usr/share/doc/asf-mapready/changelog
	echo " -- Jenkins <uso@uaf.alaska.edu>  $$(date +'%a, %d %b %Y %T %z')" \
		>> asf-mapready/usr/share/doc/asf-mapready/changelog
	gzip -9 asf-mapready/usr/share/doc/asf-mapready/changelog
	# I'm glad we got that over with.
	cp debian/copyright asf-mapready/usr/share/doc/asf-mapready/copyright
	strip asf-mapready//usr/local/bin/*
	sed -e "/Version:/s/$$/ $(DEBIAN_VERSION)/" \
	    -e "/Installed-Size:/s|$$| $$(expr $$(expr $$(du -bs asf-mapready/usr/local|cut -f1) + 512) / 1024)|" \
		-e "s/\(Architecture:\)/\1 $$(dpkg --print-architecture)/" \
		<debian/control \
		>asf-mapready/DEBIAN/control
	fakeroot dpkg-deb --build asf-mapready
ifneq (/usr/local,/usr)
	@ echo ""
	@ echo -n "          "
	@ echo "============================================================"
	@ echo -n "          "
	@ echo "Be advised that debian packages should always reside in"
	@ echo -n "          "
	@ echo "/usr/bin, so it is recommended to run ./configure"
	@ echo -n "          "
	@ echo "--prefix=/usr before again running make deb."
	@ echo -n "          "
	@ echo "============================================================"
	@ echo ""
endif
	@ echo ""
	@ echo "  XXXXXXXXXXXX Debian Package Created! XXXXXXXXXXX"
	@ echo ""

# make everything!
all: mapready oldtools req ap stp c2v jpl create_thumbs data_qc

ifneq ($(DESTDIR),"")
	DESTDIR += "/"
endif
BINDIR = $(DESTDIR)/usr/local/bin
LIBDIR = $(DESTDIR)/usr/local/lib
SHAREDIR = $(DESTDIR)/usr/local/share/asf_tools

# This should install any & all package(s) that have been built
install:
	mkdir -p "$(BINDIR)" "$(LIBDIR)" "$(SHAREDIR)"
	-cp -R $(S_BINDIR)/* "$(BINDIR)"
	-cp -R $(S_LIBDIR)/* "$(LIBDIR)"
	-cp -R $(S_SHAREDIR)/* "$(SHAREDIR)"
	chmod 755 "$(BINDIR)"/*
	if [ -d "$(SHAREDIR)/projections" ]; then \
		chmod a+x "$(SHAREDIR)"/projections ; \
		chmod a+x "$(SHAREDIR)"/projections/* ; \
	fi
	if [ -d "$(SHAREDIR)/asf_mapready" ]; then \
		chmod a+x "$(SHAREDIR)"/asf_mapready ; \
	fi
	if [ -d "$(SHAREDIR)/proj" ]; then \
		chmod a+x "$(SHAREDIR)"/proj ; \
	fi
	-cp -R $(S_DOCDIR)/*.pdf "$(SHAREDIR)/doc"
	-cp -R $(S_DOCDIR)/*.txt "$(SHAREDIR)/doc"
	chmod -R a+r "$(SHAREDIR)"
	chmod -R ug+w "$(SHAREDIR)"

# Prep the directory tree for the build
mkdirs_for_build:
	mkdir -p $(S_LIBDIR)
	mkdir -p $(S_BINDIR)
	mkdir -p $(S_DOCDIR)
	mkdir -p man/cat1
	mkdir -p man/man1

clean:
	rm -rf $(S_LIBDIR)
	rm -rf $(S_BINDIR)
	rm -rf $(S_DOCDIR)
	rm -rf $(S_SHAREDIR)
	rm -rf share
	rm -rf man
	rm -f Makefile.old
	-rm mapready*.tar.gz
	-rm -rf $(IGNORED_INCLUDES)
	-rm -rf asf-mapready
	-rm asf-mapready.deb
	$(foreach MODULE, $(ALL_MODULES), \
		$(MAKE) clean -C $(MODULE) &&) true

# Internal prep for release builds
release:
	cd make_support; $(MAKE); ./makemake linux; cd ..
