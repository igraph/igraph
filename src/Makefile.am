# This is to make sure that the headers get built before they are included
BUILT_SOURCES = foreign-ncol-parser.h foreign-lgl-parser.h

# This is needed to ensure that yacc (or bison) builds the header files
# Unfortunately this is not the default behaviour in MinGW/MSYS
AM_YFLAGS = -d

lib_LTLIBRARIES         = libigraph.la

libigraph_la_SOURCES       = basic_query.c games.c cocitation.c iterators.c \
			     structural_properties.c components.c layout.c \
			     structure_generators.c conversion.c \
			     measure_dynamics.c type_indexededgelist.c \
			     error.c interrupt.c other.c foreign.c random.c \
			     attributes.c \
			     memory.h random.h \
			     foreign-ncol-parser.y foreign-ncol-lexer.l \
			     foreign-lgl-parser.y foreign-lgl-lexer.l \
			     foreign-pajek-parser.y foreign-pajek-lexer.l \
			     foreign-gml-parser.y foreign-gml-lexer.l \
			     dqueue.c heap.c igraph_heap.c igraph_stack.c \
			     igraph_strvector.c igraph_trie.c matrix.c \
			     vector.c vector_ptr.c memory.c adjlist.c \
			     visitors.c igraph_grid.c atlas.c topology.c \
			     motifs.c progress.c atlas-edges.h operators.c \
			     igraph_psumtree.c array.c igraph_hashtable.c \
			     foreign-graphml.c flow.c igraph_buckets.c \
			     NetDataTypes.cpp NetRoutines.cpp clustertool.cpp \
			     pottsmodel_2.cpp spectral_properties.c cores.c \
			     igraph_set.c cliques.c revolver_cit.c revolver_grow.c\
			     walktrap.cpp walktrap_heap.cpp \
			     walktrap_graph.cpp walktrap_communities.cpp \
			     spmatrix.c community.c fast_community.c \
			     evolver_cit.c gml_tree.c cattributes.c \
			     revolver_ml_cit.c

libigraph_la_CFLAGS        = -fno-inline -Wall -g
libigraph_la_CXXFLAGS	   = -fno-inline -Wall -g
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm $(XML2_LIBS)

includedir = $(prefix)/include/igraph
include_HEADERS = igraph.h types.h attributes.h error.h interrupt.h \
		  NetDataTypes.h NetRoutines.h pottsmodel_2.h \
		  igraph_pmt.h igraph_pmt_off.h vector.h gml_tree.h \
		  vector.pmt walktrap_graph.h walktrap_communities.h \
		  walktrap_heap.h matrix.h matrix.pmt stack.h stack.pmt \
		  dqueue.pmt dqueue.h heap.pmt heap.h

MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h
