# This is to make sure that the headers get built before they are included
BUILT_SOURCES = foreign-ncol-parser.h foreign-lgl-parser.h

# This is needed to ensure that yacc (or bison) builds the header files
# Unfortunately this is not the default behaviour in MinGW/MSYS
AM_YFLAGS = -d

lib_LTLIBRARIES         = libigraph.la

if INTERNAL_ARPACK
  ARPACK = arpack/dgetv0.f arpack/dmout.f arpack/dsaitr.f arpack/dsapps.f \
         arpack/dsaup2.f arpack/dsaupd.f arpack/dsconv.f arpack/dseigt.f \
         arpack/dsesrt.f arpack/dseupd.f arpack/dsgets.f arpack/dsortr.f \
         arpack/dstats.f arpack/dstqrb.f arpack/dvout.f arpack/ivout.f \
	 arpack/dnaitr.f arpack/dnapps.f arpack/dnaup2.f arpack/dnaupd.f \
	 arpack/dnconv.f arpack/dsortc.f arpack/dneigh.f arpack/dngets.f \
	 arpack/dstatn.f arpack/dlaqrb.f arpack/dneupd.f \
         arpack/second.f \
	 arpack/dgeqr2.f  arpack/dlarnv.f  arpack/dlassq.f \
	 arpack/dlacpy.f  arpack/dlartg.f  arpack/dorm2r.f \
	 arpack/dlae2.f   arpack/dlaruv.f  arpack/dsteqr.f \
	 arpack/dlaev2.f  arpack/dlascl.f  arpack/dlanst.f \
	 arpack/dlaset.f  arpack/xerbla.f \
	 arpack/dlarf.f   arpack/dlasr.f   arpack/dlapy2.f \
	 arpack/dlarfg.f  arpack/dlasrt.f  arpack/dlabad.f \
	 arpack/dlanhs.f  arpack/dtrevc.f  arpack/dlanv2.f \
         arpack/dlaln2.f  arpack/dladiv.f  arpack/dtrsen.f \
         arpack/dlahqr.f  arpack/dlacon.f  arpack/dtrsyl.f \
         arpack/dtrexc.f  arpack/dlange.f  arpack/dlaexc.f \
         arpack/dlasy2.f  arpack/dlarfx.f \
	 arpack/daxpy.f  arpack/ddot.f   arpack/dger.f   \
	 arpack/dscal.f  arpack/lsame.f  arpack/dcopy.f  \
	 arpack/dgemv.f  arpack/dnrm2.f  arpack/dswap.f  \
	 arpack/drot.f   arpack/idamax.f arpack/dtrmm.f  \
	 arpack/dasum.f
  noinst_LTLIBRARIES  = libdlamch.la
  libdlamch_la_SOURCES = arpack/dlamch.f
  libdlamch_la_FFLAGS = -ffloat-store
  ARPACKLAPACK_DLAMCH = libdlamch.la
else
  ARPACK = arpack/xerbla.f
  ARPACKLAPACK_DLAMCH =
endif

if INTERNAL_LAPACK
  LAPACK = lapack/ldgemm.f lapack/llsame.f lapack/lxerbla.f
else
  LAPACK = 
endif

libigraph_la_SOURCES       = scg.c \
			     basic_query.c games.c cocitation.c iterators.c \
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
			     forestfire.c $(ARPACK) $(LAPACK) \
			     arpack.c centrality.c drl_layout.cpp drl_parse.cpp \
                             drl_graph.cpp DensityGrid.cpp \
			     gengraph_box_list.cpp gengraph_degree_sequence.cpp \
			     gengraph_graph_molloy_hash.cpp \
			     gengraph_graph_molloy_optimized.cpp \
			     gengraph_mr-connected.cpp gengraph_powerlaw.cpp \
			     gengraph_random.cpp decomposition.c bipartite.c \
			     scg_approximate_methods.c scg_exact_scg.c \
			     scg_grouping.c scg_kmeans.c scg_optimal_method.c \
			     scg_utils.c scg_headers.h igraph_complex.c

libigraph_la_CFLAGS        = -Wall
libigraph_la_CXXFLAGS	   = -Wall
libigraph_la_FFLAGS        = -ffixed-line-length-132
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm $(XML2_LIBS) $(FLIBS) $(ARPACKLAPACK_DLAMCH)

includedir = $(prefix)/include/igraph
include_HEADERS = igraph.h types.h attributes.h error.h interrupt.h \
		  NetDataTypes.h NetRoutines.h pottsmodel_2.h \
		  igraph_pmt.h igraph_pmt_off.h vector.h gml_tree.h \
		  vector.pmt walktrap_graph.h walktrap_communities.h \
		  walktrap_heap.h matrix.h matrix.pmt stack.h stack.pmt \
		  dqueue.pmt dqueue.h heap.pmt heap.h arpack.h array.pmt \
		  array.h arpack_internal.h f2c.h igraph_math.h igraph_complex.h


MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h
