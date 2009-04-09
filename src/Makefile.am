# This is to make sure that the headers get built before they are included
BUILT_SOURCES = foreign-ncol-parser.h foreign-lgl-parser.h

# This is needed to ensure that yacc (or bison) builds the header files
# Unfortunately this is not the default behaviour in MinGW/MSYS
AM_YFLAGS = -d

lib_LTLIBRARIES         = libigraph.la

if INTERNAL_LAPACK
  LAPACK = lapack/dgeqr2.c  lapack/dlarnv.c  lapack/dlassq.c \
	lapack/dlacpy.c  lapack/dlartg.c  lapack/dorm2r.c \
	lapack/dlae2.c   lapack/dlaruv.c  lapack/dsteqr.c \
	lapack/dlaev2.c  lapack/dlascl.c  lapack/dlanst.c \
	lapack/dlaset.c  lapack/xerbla.c \
	lapack/dlarf.c   lapack/dlasr.c   lapack/dlapy2.c \
	lapack/dlarfg.c  lapack/dlasrt.c  lapack/dlabad.c \
	lapack/dlanhs.c  lapack/dtrevc.c  lapack/dlanv2.c \
        lapack/dlaln2.c  lapack/dladiv.c  lapack/dtrsen.c \
        lapack/dlahqr.c  lapack/dlacon.c  lapack/dtrsyl.c \
        lapack/dtrexc.c  lapack/dlange.c  lapack/dlaexc.c \
        lapack/dlasy2.c  lapack/dlarfx.c 
  noinst_LTLIBRARIES  = libdlamch.la
  libdlamch_la_SOURCES = lapack/dlamch.c
  libdlamch_la_CFLAGS = $(FLOATSTORE)
  LAPACK_DLAMCH = libdlamch.la
else 
  LAPACK = lapack/xerbla.c
  LAPACK_DLAMCH = 
endif

if INTERNAL_ARPACK
  ARPACK = arpack/dgetv0.c arpack/dmout.c arpack/dsaitr.c arpack/dsapps.c \
         arpack/dsaup2.c arpack/dsaupd.c arpack/dsconv.c arpack/dseigt.c \
         arpack/dsesrt.c arpack/dseupd.c arpack/dsgets.c arpack/dsortr.c \
         arpack/dstats.c arpack/dstqrb.c arpack/dvout.c arpack/ivout.c \
	 arpack/dnaitr.c arpack/dnapps.c arpack/dnaup2.c arpack/dnaupd.c \
	 arpack/dnconv.c arpack/dsortc.c arpack/dneigh.c arpack/dngets.c \
	 arpack/dstatn.c arpack/dlaqrb.c arpack/dneupd.c \
         arpack/second.c
else
  ARPACK =
endif

if INTERNAL_BLAS
  BLAS = blas/daxpy.c  blas/ddot.c   blas/dger.c   blas/dscal.c  blas/lsame.c \
	blas/dcopy.c  blas/dgemv.c  blas/dnrm2.c  blas/dswap.c   blas/drot.c  \
        blas/idamax.c blas/dtrmm.c  blas/dasum.c
else
  BLAS =
endif

F2C = f2c/d_sign.c f2c/etime_.c  f2c/pow_dd.c f2c/pow_di.c f2c/s_cmp.c \
	f2c/s_copy.c f2c/d_lg10.c f2c/i_dnnt.c

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
			     forestfire.c $(ARPACK) $(BLAS) $(LAPACK) $(F2C) \
			     arpack.c centrality.c drl_layout.cpp drl_parse.cpp \
                             drl_graph.cpp DensityGrid.cpp \
			     gengraph_box_list.cpp gengraph_degree_sequence.cpp \
			     gengraph_graph_molloy_hash.cpp \
			     gengraph_graph_molloy_optimized.cpp \
			     gengraph_mr-connected.cpp gengraph_powerlaw.cpp \
			     gengraph_random.cpp bipartite.c \
			     drl_layout_3d.cpp drl_graph_3d.cpp \
			     DensityGrid_3d.cpp

libigraph_la_CFLAGS        =
libigraph_la_CXXFLAGS	   =
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm $(XML2_LIBS) $(FLIBS) $(LAPACK_DLAMCH)

includedir = $(prefix)/include/igraph
include_HEADERS = igraph.h types.h attributes.h error.h interrupt.h \
		  NetDataTypes.h NetRoutines.h pottsmodel_2.h \
		  igraph_pmt.h igraph_pmt_off.h vector.h gml_tree.h \
		  vector.pmt walktrap_graph.h walktrap_communities.h \
		  walktrap_heap.h matrix.h matrix.pmt stack.h stack.pmt \
		  dqueue.pmt dqueue.h heap.pmt heap.h arpack.h array.pmt \
		  array.h arpack_internal.h f2c.h drl_layout.h \
		  drl_parse.h drl_graph.h DensityGrid.h drl_Node.h \
		  gengraph_box_list.h gengraph_header.h \
		  gengraph_definitions.h gengraph_powerlaw.h \
		  gengraph_degree_sequence.h gengraph_qsort.h \
		  gengraph_graph_molloy_hash.h gengraph_random.h \
		  gengraph_graph_molloy_optimized.h gengraph_vertex_cover.h \
		  gengraph_hash.h igraph_math.h \
		  drl_graph_3d.h drl_layout_3d.h drl_Node_3d.h \
		  DensityGrid_3d.h


MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h
