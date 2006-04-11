pyexec_LTLIBRARIES = igraph.la

igraph_la_SOURCES       = igraphmodule.c error.c error.h convert.c convert.h \
                          vertexobject.c vertexobject.h \
                          vertexseqobject.c vertexseqobject.h \
                          edgeobject.c edgeobject.h \
			  edgeseqobject.c edgeseqobject.h \
			  graphobject.c graphobject.h \
			  common.c common.h
igraph_la_CFLAGS        = $(PYTHON_INCLUDES) -I$(top_srcdir)/include
igraph_la_LDFLAGS       = -module $(PYTHON_LDFLAGS) -L../../src
igraph_la_LIBADD        = -lm $(PYTHON_LIBS) -ligraph

EXTRA_DIST = test.py test.sh

## We don't need libtool's versioning information in the Python
## library, so remove them
install-data-hook:
	rm -f $(pyexecdir)/igraph.so $(pyexecdir)/igraph.so.0 $(pyexecdir)/igraph.so.0.0
	rm -f $(pyexecdir)/igraph.a $(pyexecdir)/igraph.la
	mv $(pyexecdir)/igraph.so.0.0.0 $(pyexecdir)/igraph.so

MAINTAINERCLEANFILES = Makefile.in
