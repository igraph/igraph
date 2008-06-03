bin_PROGRAMS = igraph

STIMULUS = $(top_srcdir)/tools/stimulus.py
FUNCTIONS_DEF = $(top_srcdir)/interfaces/functions.def

interface.c: interface.c.in types.def $(STIMULUS) $(FUNCTIONS_DEF)
	$(STIMULUS) -f $(FUNCTIONS_DEF) \
	  -i interface.c.in \
	  -o interface.c \
	  -t types.def \
          -l Shell

funcs.txt: interface.c.in $(STIMULUS) $(FUNCTIONS_DEF)
	$(STIMULUS) -f $(FUNCTIONS_DEF) -o funcs.txt -l ShellLn

igraph_CPPFLAGS = -I $(top_srcdir)/include/
igraph_CFLAGS = -Wall -O2
igraph_LDFLAGS = -L$(top_srcdir)/src/.libs -ligraph
igraph_SOURCES = interface.c

install-exec-hook: funcs.txt
	for i in `cat funcs.txt`; do ln -sf igraph $(bindir)/$$i; done
