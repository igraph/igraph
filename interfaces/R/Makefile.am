
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
	operators.c interrupt.c igraph_psumtree.c array.c \
	igraph_hashtable.c foreign-graphml.c flow.c igraph_buckets.c \
        spectral_properties.c cores.c cliques.c igraph_set.c \
	evolver.c community.c spmatrix.c

RSNAMES_F = 

RSNAMES_CPP = NetDataTypes.cpp clustertool.cpp NetRoutines.cpp \
	pottsmodel_2.cpp walktrap.cpp walktrap_heap.cpp walktrap_graph.cpp \
	walktrap_communities.cpp fastcommunity_mh.cc

RSNAMES_H = memory.h random.h atlas-edges.h NetDataTypes.h \
	NetRoutines.h pottsmodel_2.h walktrap_communities.h \
	walktrap_heap.h walktrap_graph.h fastcommunity_vektor.h \
	fastcommunity_maxheap.h

RSNAMES_Y = foreign-ncol-parser.y foreign-lgl-parser.y foreign-pajek-parser.y

RSNAMES_L = foreign-ncol-lexer.l foreign-lgl-lexer.l foreign-pajek-lexer.l

RSCNAMES = $(RSNAMES_C) $(RSNAMES_H) $(RSNAMES_CPP) $(RSNAMES_F)

RSCNAMES2 = $(RSNAMES_Y:.y=.c) $(RSNAMES_Y:.y=.h) \
	$(RSNAMES_L:.l=.c)

RSCNAMES_ORIGPLACE = $(patsubst %,$(top_srcdir)/src/%,$(RSCNAMES)) \
	$(patsubst %,$(top_builddir)/src/%,$(RSCNAMES2))

RTOPHEADERS = $(top_srcdir)/include/attributes.h \
	$(top_srcdir)/include/error.h \
	$(top_srcdir)/include/igraph.h \
	$(top_srcdir)/include/types.h \
	$(top_srcdir)/include/interrupt.h

RSRCSOURCES = $(top_srcdir)/interfaces/R/src/rinterface.c \
	$(top_srcdir)/interfaces/R/src/Makevars.in \
	$(top_srcdir)/interfaces/R/src/Makevars.win \
	$(top_srcdir)/interfaces/R/src/config.h.in  \
	$(top_srcdir)/interfaces/R/src/config.h

# ../../interfaces/R/src/rinterface.c: config.def ../functions.def ../../interfaces/R/src/rinterface.c.in
# 	../../tools/inger.py -c config.def -f ../functions.def

# igraph/R/auto.R: config.def ../functions.def igraph/R/auto.R.in
# 	../../tools/inger.py -c config.def -f ../functions.def

# igraph/NAMESPACE: config.def ../functions.def igraph/NAMESPACE.in
# 	../../tools/inger.py -c config.def -f ../functions.def

igraph_$(VERSION).tar.gz: $(RSCNAMES_ORIGPLACE) $(RSRCSOURCES) $(RRFILES)
	if [ "x$(top_srcdir)" != "x$(top_builddir)" ]; then cp -r $(top_srcdir)/interfaces/R/igraph .; fi
	cp $(RSCNAMES_ORIGPLACE) igraph/src/
	cp $(RTOPHEADERS) igraph/src
	cp $(RSRCSOURCES) igraph/src
	cp DESCRIPTION igraph/
	R CMD build igraph

CLEANFILES = igraph/src/*

