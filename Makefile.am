
if DOCS
   MAYBE_DOCS = doc
endif

SUBDIRS=src interfaces debian tests $(MAYBE_DOCS)
DOC_FILES = $(top_builddir)/doc/html \
	$(top_builddir)/doc/igraph-docs.info \
	$(top_builddir)/doc/igraph-docs.pdf

pkgconfigdir = @libdir@/pkgconfig
pkgconfig_DATA = igraph.pc

EXTRA_DIST = igraph.pc $(top_srcdir)/include/* $(DOC_FILES) \
	examples 

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
