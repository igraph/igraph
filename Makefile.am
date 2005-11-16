SUBDIRS=src interfaces debian tests

pkgconfigdir = @libdir@/pkgconfig
pkgconfig_DATA = igraph.pc

EXTRA_DIST = igraph.pc include/*

MAINTAINERCLEANFILES = Makefile.in

## to make sure make doc will always call doxygen rule
## and deb will generate Debian packages
.PHONY: doc deb debcvs

doc:
	doxygen Doxyfile

deb: clean
	debian/prepare
	fakeroot debian/rules binary

debcvs: clean
	debian/prepare cvs
	fakeroot debian/rules binary
	
test: all
	for i in interfaces/*; do \
		if [ -d $$i ]; then \
		  $$i/test.sh; \
	 	fi; \
	done

CLEANFILES=doc/*

