SUBDIRS=src interfaces debian

pkgconfigdir = @libdir@/pkgconfig
pkgconfig_DATA = igraph.pc

EXTRA_DIST = igraph.pc include/*

doc: 
	doxygen Doxyfile

deb:
	debian/prepare
	fakeroot debian/rules binary

debcvs:
	debian/prepare cvs
	fakeroot debian/rules binary
	
test: all
	for i in interfaces/*; do \
		if [ -d $$i ]; then \
		  $$i/test.sh; \
	 	fi; \
	done

CLEANFILES=doc/*

## to make sure make doc will always call doxygen rule
## and deb will generate Debian packages
.PHONY: doc deb
