bin_PROGRAMS = igraph

interface.c: interface.c.in ../functions.def types.def ../../tools/stimulus.py
	../../tools/stimulus.py \
	  -f ../functions.def \
	  -i interface.c.in \
	  -o interface.c \
	  -t types.def \
          -l Shell

igraph_CPPFLAGS = -I $(top_srcdir)/include/
igraph_LDFLAGS = -L$(top_srcdir)/src/.libs -ligraph
igraph_SOURCES = interface.c
