# This is to make sure that the headers get built before they are included
BUILT_SOURCES = foreign-ncol-parser.h foreign-lgl-parser.h

# This is needed to ensure that yacc (or bison) builds the header files
# Unfortunately this is not the default behaviour in MinGW/MSYS
AM_YFLAGS = -d

if INTERNAL_LAPACK
  LAPACK = lapack/lapack_igraph.f lapack/print.c
else 
  LAPACK = lapack/print.c
endif

if INTERNAL_ARPACK
  ARPACK = arpack/dgetv0.f arpack/dmout.f arpack/dsaitr.f arpack/dsapps.f \
         arpack/dsaup2.f arpack/dsaupd.f arpack/dsconv.f arpack/dseigt.f \
         arpack/dsesrt.f arpack/dseupd.f arpack/dsgets.f arpack/dsortr.f \
         arpack/dstats.f arpack/dstqrb.f arpack/dvout.f arpack/ivout.f \
         arpack/second.f
else
  ARPACK =
endif

if INTERNAL_BLAS
  BLAS = blas/blas_igraph.f
else
  BLAS =
endif

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
			     evolver_cit.c gml_tree.c \
			     bliss_bignum.hh bliss_kstack.hh \
			     bliss_defs.hh bliss_orbit.cc \
			     bliss_eqrefhash.cc bliss_orbit.hh \
			     bliss_eqrefhash.hh bliss_partition.cc \
			     bliss_graph.cc bliss_partition.hh \
			     bliss_graph.hh bliss_timer.cc \
			     bliss_heap.cc bliss_timer.hh \
			     bliss_heap.hh bliss_utils.cc \
			     bliss_kqueue.hh bliss_utils.hh bliss.cc \
			     cattributes.c revolver_ml_cit.c zeroin.c bfgs.c math.c \
			     forestfire.c $(ARPACK) $(BLAS) $(LAPACK) \
			     arpack.c

libigraph_la_CFLAGS        = 
libigraph_la_CXXFLAGS	   = 
libigraph_la_FFLAGS        =
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm $(XML2_LIBS) $(FLIBS)

includedir = $(prefix)/include/igraph
include_HEADERS = igraph.h types.h attributes.h error.h interrupt.h \
		  NetDataTypes.h NetRoutines.h pottsmodel_2.h \
		  igraph_pmt.h igraph_pmt_off.h vector.h gml_tree.h \
		  vector.pmt walktrap_graph.h walktrap_communities.h \
		  walktrap_heap.h matrix.h matrix.pmt stack.h stack.pmt \
		  dqueue.pmt dqueue.h heap.pmt heap.h

MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h
