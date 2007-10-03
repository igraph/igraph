bin_PROGRAMS = igraph

interface.c: interface.c.in ../functions.def types.def ../../tools/stimulus.py
	../../tools/stimulus.py \
	  -f ../functions.def \
	  -i interface.c.in \
	  -o interface.c \
	  -t types.def \
          -l Shell

funcs.txt: interface.c.in ../functions.def ../../tools/stimulus.py
	../../tools/stimulus.py \
	-f ../functions.def \
	-o funcs.txt \
	-l ShellLn

igraph_CPPFLAGS = -I $(top_srcdir)/include/
igraph_CFLAGS = -Wall -O2
igraph_LDFLAGS = -L$(top_srcdir)/src/.libs -ligraph
igraph_SOURCES = interface.c

install-exec-hook: funcs.txt
	for i in `cat funcs.txt`; do ln -sf igraph $(bindir)/$$i; done
