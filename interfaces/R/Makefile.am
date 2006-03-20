
all: igraph_$(VERSION).tar.gz

RRFILES = $(wildcard igraph/R/*.R) \
	  $(wildcard igraph/man/*.Rd) \
	  igraph/NAMESPACE igraph/DESCRIPTION

RSNAMES_C = basic_query.c games.c cocitation.c iterators.c \
	structural_properties.c components.c layout.c \
	structure_generators.c conversion.c \
	measure_dynamics.c type_indexededgelist.c \
	error.c other.c foreign.c random.c \
	attributes.c dqueue.c heap.c igraph_stack.c \
	igraph_strvector.c igraph_trie.c matrix.c \
	vector.c vector_ptr.c memory.c adjlist.c visitors.c \
	igraph_grid.c atlas.c topology.c motifs.c progress.c \
	operators.c interrupt.c

RSNAMES_H = memory.h random.h atlas-edges.h

RSNAMES_Y = foreign-ncol-parser.y foreign-lgl-parser.y foreign-pajek-parser.y

RSNAMES_L = foreign-ncol-lexer.l foreign-lgl-lexer.l foreign-pajek-lexer.l

RSCNAMES = $(RSNAMES_C) $(RSNAMES_H) \
	$(RSNAMES_Y:.y=.c) $(RSNAMES_Y:.y=.h) \
	$(RSNAMES_L:.l=.c)

RSCNAMES_ORIGPLACE = $(patsubst %,$(top_srcdir)/src/%,$(RSCNAMES))

RTOPHEADERS = $(top_srcdir)/include/attributes.h \
	$(top_srcdir)/include/error.h \
	$(top_srcdir)/include/igraph.h \
	$(top_srcdir)/include/types.h \
	$(top_srcdir)/include/interrupt.h

RSRCSOURCES = src/rinterface.c src/Makevars src/Makevars.win \
	src/config.h.in src/config.h

igraph_$(VERSION).tar.gz: $(RSCNAMES_ORIGPLACE) $(RSRCSOURCES) $(RRFILES)
	cp $(RSCNAMES_ORIGPLACE) igraph/src/
	cp $(RTOPHEADERS) igraph/src
	cp $(RSRCSOURCES) igraph/src
	cp DESCRIPTION igraph/
	R CMD build igraph

CLEANFILES = igraph/src/*

