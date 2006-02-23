
if DOCS
   MAYBE_DOCS = doc
endif

SUBDIRS=src interfaces debian tests $(MAYBE_DOCS)
DOC_FILES = doc/html doc/igraph-docs.info doc/igraph-docs.pdf

pkgconfigdir = @libdir@/pkgconfig
pkgconfig_DATA = igraph.pc

EXTRA_DIST = igraph.pc include/* $(DOC_FILES) \
	examples tests

MAINTAINERCLEANFILES = Makefile.in

## to make sure make doc will always call doxygen rule
## and deb will generate Debian packages
.PHONY: dox deb debcvs

dox:
	doxygen Doxyfile

deb: clean
	debian/prepare
	fakeroot debian/rules binary

debcvs: clean
	debian/prepare cvs
	fakeroot debian/rules binary

test: all
	for i in interfaces/*; do \
		if [ -x $$i/test.sh ]; then \
		  $$i/test.sh; \
	 	fi; \
	done

CLEANFILES=
