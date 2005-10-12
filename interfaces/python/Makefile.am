lib_LTLIBRARIES         = igraph.la

igraph_la_SOURCES       = igraphmodule.c
igraph_la_CFLAGS        = $(PYTHON_INCLUDES) -I$(top_srcdir)/include
igraph_la_LDFLAGS       = -module
igraph_la_LIBADD        = -lm $(PYTHON_LIBS) -ligraph

## need to move the generated Python libraries to their proper place.
## unfortunately Automake does not have any native support to define
## destination directories for libraries
install-exec-hook:
	mkdir --parents $(DESTDIR)$(pythondir)
	mv `readlink -f $(DESTDIR)$(libdir)/igraph.so` $(DESTDIR)$(pythondir)/igraph.so
	rm $(DESTDIR)$(libdir)/igraph.so.* $(DESTDIR)$(libdir)/igraph.a $(DESTDIR)$(libdir)/igraph.la
