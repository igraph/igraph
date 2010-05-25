# This is to make sure that the headers get built before they are included
BUILT_SOURCES = foreign-ncol-parser.h foreign-lgl-parser.h \
		foreign-dl-parser.h foreign-gml-parser.h \
		foreign-pajek-parser.h

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
  libdlamch_la_CFLAGS = $(FLOATSTORE) -I$(top_srcdir)/include
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

CS = cs/cs_add.c       cs/cs_happly.c    cs/cs_pvec.c         \
	cs/cs_amd.c       cs/cs_house.c     cs/cs_qr.c        \
	cs/cs_chol.c      cs/cs_ipvec.c     cs/cs_qrsol.c     \
	cs/cs_cholsol.c   cs/cs_leaf.c      cs/cs_randperm.c  \
	cs/cs_compress.c  cs/cs_load.c      cs/cs_reach.c     \
	cs/cs_counts.c    cs/cs_lsolve.c    cs/cs_scatter.c   \
	cs/cs_cumsum.c    cs/cs_ltsolve.c   cs/cs_scc.c       \
	cs/cs_dfs.c       cs/cs_lu.c        cs/cs_schol.c     \
	cs/cs_dmperm.c    cs/cs_lusol.c     cs/cs_spsolve.c   \
	cs/cs_droptol.c   cs/cs_malloc.c    cs/cs_sqr.c       \
	cs/cs_dropzeros.c cs/cs_maxtrans.c  cs/cs_symperm.c   \
	cs/cs_dupl.c      cs/cs_multiply.c  cs/cs_tdfs.c      \
	cs/cs_entry.c     cs/cs_norm.c      cs/cs_transpose.c \
	cs/cs_ereach.c    cs/cs_permute.c   cs/cs_updown.c    \
	cs/cs_etree.c     cs/cs_pinv.c      cs/cs_usolve.c    \
	cs/cs_fkeep.c     cs/cs_post.c      cs/cs_util.c      \
	cs/cs_gaxpy.c     cs/cs_print.c     cs/cs_utsolve.c   \
	cs/cs.h	    cs/UFconfig.h

HEADERS_PRIVATE = atlas-edges.h \
		bliss_bignum.hh		bliss_defs.hh \
		bliss_eqrefhash.hh	bliss_graph.hh \
		bliss_heap.hh		bliss_kqueue.hh \
		bliss_kstack.hh		bliss_orbit.hh \
		bliss_partition.hh	bliss_timer.hh \
		bliss_utils.hh \
		NetDataTypes.h		NetRoutines.h \
		pottsmodel_2.h \
		igraph_gml_tree.h \
		walktrap_graph.h	walktrap_communities.h \
		walktrap_heap.h \
		igraph_blas_internal.h \
		igraph_arpack_internal.h	igraph_f2c.h \
		igraph_math.h \
		drl_layout.h drl_parse.h drl_graph.h \
		drl_graph_3d.h		drl_layout_3d.h \
		drl_Node.h	drl_Node_3d.h \
		DensityGrid.h	DensityGrid_3d.h \
		vector.pmt matrix.pmt stack.pmt dqueue.pmt heap.pmt array.pmt \
		igraph_types_internal.h \
		foreign-dl-header.h bignum.h bigint.h \
		gengraph_box_list.h gengraph_definitions.h \
		gengraph_degree_sequence.h gengraph_graph_molloy_hash.h \
		gengraph_graph_molloy_optimized.h \
		gengraph_hash.h gengraph_header.h gengraph_powerlaw.h \
		gengraph_qsort.h gengraph_random.h gengraph_vertex_cover.h

HEADERS_PUBLIC =../include/igraph.h 		../include/igraph_memory.h    \
		../include/igraph_random.h 	../include/igraph_types.h     \
		../include/igraph_vector.h 	../include/igraph_matrix.h    \
		../include/igraph_array.h 	../include/igraph_dqueue.h    \
		../include/igraph_stack.h 	../include/igraph_heap.h      \
		../include/igraph_arpack.h 	../include/igraph_interrupt.h \
		../include/igraph_attributes.h 	../include/igraph_error.h     \
		../include/igraph_pmt.h 	../include/igraph_pmt_off.h   \
		../include/igraph_adjlist.h	../include/igraph_iterators.h \
		../include/igraph_bipartite.h	../include/igraph_layout.h    \
		../include/igraph_centrality.h	../include/igraph_motifs.h    \
		../include/igraph_cliques.h	../include/igraph_neighborhood.h \
		../include/igraph_cocitation.h	../include/igraph_nongraph.h  \
		../include/igraph_community.h	../include/igraph_operators.h \
		../include/igraph_components.h	../include/igraph_paths.h     \
		../include/igraph_constructors.h ../include/igraph_progress.h \
		../include/igraph_conversion.h	../include/igraph_revolver.h  \
		../include/igraph_datatype.h	../include/igraph_structural.h\
		../include/igraph_flow.h	../include/igraph_topology.h  \
		../include/igraph_foreign.h	../include/igraph_transitivity.h \
		../include/igraph_games.h	../include/igraph_visitor.h   \
		../include/igraph_interface.h	../include/igraph_adjlist.h   \
		../include/igraph_progress.h	../include/igraph_constants.h \
		../include/igraph_vector_pmt.h	../include/igraph_matrix_pmt.h\
		../include/igraph_array_pmt.h	../include/igraph_dqueue_pmt.h\
		../include/igraph_stack_pmt.h	../include/igraph_heap_pmt.h  \
		../include/igraph_vector_ptr.h	../include/igraph_spmatrix.h  \
		../include/igraph_strvector.h	../include/igraph_psumtree.h  \
		../include/igraph_sparsemat.h   ../include/igraph_mixing.h    \
		../include/igraph_version.h     ../include/igraph_blas.h


libigraph_la_SOURCES       = $(HEADERS_PRIVATE) \
			     basic_query.c games.c cocitation.c iterators.c \
			     structural_properties.c components.c layout.c \
			     structure_generators.c conversion.c \
			     measure_dynamics.c type_indexededgelist.c \
			     error.c interrupt.c other.c foreign.c random.c \
			     attributes.c \
			     foreign-ncol-parser.y foreign-ncol-lexer.l \
			     foreign-lgl-parser.y foreign-lgl-lexer.l \
			     foreign-pajek-parser.y foreign-pajek-lexer.l \
			     foreign-gml-parser.y foreign-gml-lexer.l \
			     dqueue.c heap.c igraph_heap.c igraph_stack.c \
			     igraph_strvector.c igraph_trie.c matrix.c \
			     vector.c vector_ptr.c memory.c adjlist.c \
			     visitors.c igraph_grid.c atlas.c topology.c \
			     motifs.c progress.c operators.c \
			     igraph_psumtree.c array.c igraph_hashtable.c \
			     foreign-graphml.c flow.c igraph_buckets.c \
			     NetDataTypes.cpp NetRoutines.cpp clustertool.cpp \
			     pottsmodel_2.cpp spectral_properties.c cores.c \
			     igraph_set.c cliques.c revolver_cit.c revolver_grow.c\
			     walktrap.cpp walktrap_heap.cpp \
			     walktrap_graph.cpp walktrap_communities.cpp \
			     spmatrix.c community.c fast_community.c \
			     evolver_cit.c gml_tree.c \
			     bliss_orbit.cc bliss_eqrefhash.cc \
			     bliss_partition.cc bliss_graph.cc bliss_timer.cc \
			     bliss_heap.cc bliss_utils.cc bliss.cc \
			     cattributes.c revolver_ml_cit.c zeroin.c bfgs.c math.c \
			     forestfire.c $(ARPACK) $(BLAS) $(LAPACK) $(F2C) \
			     blas.c arpack.c centrality.c drl_layout.cpp drl_parse.cpp \
			     drl_graph.cpp DensityGrid.cpp \
			     gengraph_box_list.cpp gengraph_degree_sequence.cpp \
			     gengraph_graph_molloy_hash.cpp \
			     gengraph_graph_molloy_optimized.cpp \
			     gengraph_mr-connected.cpp gengraph_powerlaw.cpp \
			     gengraph_random.cpp decomposition.c bipartite.c \
			     drl_layout_3d.cpp drl_graph_3d.cpp \
			     DensityGrid_3d.cpp \
			     foreign-dl-parser.y foreign-dl-lexer.l \
			     $(CS) sparsemat.c mixing.c bigint.c bignum.c \
			     version.c optimal_modularity.c \
			     igraph_fixed_vectorlist.c

libigraph_la_CFLAGS        = -I$(top_srcdir)/include -I$(top_builddir)/include -Wall
libigraph_la_CXXFLAGS	   = -I$(top_srcdir)/include -I$(top_builddir)/include -Wall
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm $(XML2_LIBS) $(FLIBS) $(LAPACK_DLAMCH)

includedir = $(prefix)/include/igraph
include_HEADERS = $(HEADERS_PUBLIC)

MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h
