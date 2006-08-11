pyexec_LTLIBRARIES = igraph.la

igraph_la_SOURCES       = igraphmodule.c error.c error.h convert.c convert.h \
                          vertexobject.c vertexobject.h \
                          vertexseqobject.c vertexseqobject.h \
                          edgeobject.c edgeobject.h \
			  edgeseqobject.c edgeseqobject.h \
			  graphobject.c graphobject.h \
			  bfsiter.c bfsiter.h \
			  common.c common.h
igraph_la_CFLAGS        = $(PYTHON_INCLUDES) -I$(top_srcdir)/include
igraph_la_LDFLAGS       = -module $(PYTHON_LDFLAGS) -L../../src
igraph_la_LIBADD        = -lm $(PYTHON_LIBS) -ligraph

EXTRA_DIST = test.py test.sh

## We don't need libtool's versioning information in the Python
## library, so remove them
install-data-hook:
	rm -f $(pyexecdir)/igraph.so $(pyexecdir)/igraph.so.0 $(pyexecdir)/igraph.so.0.0 $(pyexecdir)/igraph.0.so $(pyexecdir)/igraph.0.0.so
	rm -f $(pyexecdir)/igraph.a $(pyexecdir)/igraph.la
	if [ -f $(pyexecdir)/igraph.so.0.0.0 ]; then mv $(pyexecdir)/igraph.so.0.0.0 $(pyexecdir)/igraph.so; fi
	if [ -f $(pyexecdir)/igraph.0.0.0.so ]; then mv $(pyexecdir)/igraph.0.0.0.so $(pyexecdir)/igraph.so; fi

MAINTAINERCLEANFILES = Makefile.in
