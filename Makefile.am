

SUBDIRS=src debian tests interfaces
DOC_FILES = 	$(top_builddir)/doc/book/Makefile.in \
	$(top_builddir)/doc/igraph-docs.info \
	$(top_srcdir)/doc/Makefile.in

pkgconfigdir = @libdir@/pkgconfig
pkgconfig_DATA = igraph.pc

EXTRA_DIST = igraph.pc $(top_srcdir)/include/* $(DOC_FILES) \
	examples $(top_srcdir)/interfaces/R/Makefile.in \
	$(top_srcdir)/interfaces/R/DESCRIPTION.in

install-exec-hook:
	if test -f $(top_builddir)/src/.libs/cygigraph-0.dll ; \
		then cp $(top_builddir)/src/.libs/cygigraph-0.dll \
		$(DESTDIR)$(libdir) ; fi

MAINTAINERCLEANFILES = Makefile.in

## to make sure make doc will always call doxygen rule
## and deb will generate Debian packages
.PHONY: dox deb debcvs

dox:
	doxygen Doxyfile

deb: clean
	debian/prepare
	dpkg-buildpackage -rfakeroot

test: all
	for i in interfaces/*; do \
		if [ -x $$i/test.sh ]; then \
		  $$i/test.sh; \
	 	fi; \
	done

CLEANFILES=
