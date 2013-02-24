# This is to make sure that the headers get built before they are included
BUILT_SOURCES = foreign-ncol-parser.h foreign-lgl-parser.h \
		foreign-dl-parser.h foreign-gml-parser.h \
		foreign-pajek-parser.h

# This is needed to ensure that yacc (or bison) builds the header files
# Unfortunately this is not the default behaviour in MinGW/MSYS
AM_YFLAGS = -d

lib_LTLIBRARIES         = libigraph.la

include lapack/blas.inc
include lapack/lapack.inc
include lapack/arpack.inc
include plfit/plfit.inc

F2C =	f2c/abort_.c	  f2c/dolio.c	    f2c/r_sin.c\
	f2c/dummy.c	  f2c/dtime_.c	    f2c/iio.c	      f2c/r_sinh.c\
	f2c/backspac.c	  f2c/due.c	    f2c/ilnw.c	      f2c/r_sqrt.c\
	f2c/c_abs.c	  f2c/ef1asc_.c	    f2c/inquire.c     f2c/r_tan.c\
	f2c/c_cos.c	  f2c/ef1cmc_.c	    f2c/l_ge.c	      f2c/r_tanh.c\
	f2c/c_div.c	  f2c/endfile.c	    f2c/l_gt.c	      f2c/rdfmt.c\
	f2c/c_exp.c	  f2c/erf_.c	    f2c/l_le.c	      f2c/rewind.c\
	f2c/c_log.c	  f2c/erfc_.c	    f2c/l_lt.c	      f2c/rsfe.c\
	f2c/c_sin.c	  f2c/err.c	    f2c/lbitbits.c    f2c/rsli.c\
	f2c/c_sqrt.c	  f2c/etime_.c	    f2c/lbitshft.c    f2c/rsne.c\
	f2c/cabs.c	  f2c/exit_.c	    f2c/lread.c	      f2c/s_cat.c\
	f2c/close.c	  f2c/f77_aloc.c    f2c/lwrite.c      f2c/s_cmp.c\
	f2c/ctype.c	  f2c/f77vers.c	    f2c/s_copy.c\
	f2c/d_abs.c	  f2c/fmt.c	    f2c/open.c	      f2c/s_paus.c\
	f2c/d_acos.c	  f2c/fmtlib.c	    f2c/pow_ci.c      f2c/s_rnge.c\
	f2c/d_asin.c	  f2c/ftell_.c	    f2c/pow_dd.c      f2c/s_stop.c\
	f2c/d_atan.c	  f2c/pow_di.c      f2c/sfe.c\
	f2c/d_atn2.c	  f2c/getenv_.c	    f2c/pow_hh.c      f2c/sig_die.c\
	f2c/d_cnjg.c	  f2c/h_abs.c	    f2c/pow_ii.c      f2c/signal_.c\
	f2c/d_cos.c	  f2c/h_dim.c	    f2c/pow_ri.c      f2c/signbit.c\
	f2c/d_cosh.c	  f2c/h_dnnt.c	    f2c/pow_zi.c      f2c/sue.c\
	f2c/d_dim.c	  f2c/h_indx.c	    f2c/pow_zz.c      f2c/system_.c\
	f2c/d_exp.c	  f2c/h_len.c	    f2c/r_abs.c	      f2c/typesize.c\
	f2c/d_imag.c	  f2c/h_mod.c	    f2c/r_acos.c      f2c/uio.c\
	f2c/d_int.c	  f2c/h_nint.c	    f2c/r_asin.c      f2c/uninit.c\
	f2c/d_lg10.c	  f2c/h_sign.c	    f2c/r_atan.c      f2c/util.c\
	f2c/d_log.c	  f2c/hl_ge.c	    f2c/r_atn2.c      f2c/wref.c\
	f2c/d_mod.c	  f2c/hl_gt.c	    f2c/r_cnjg.c      f2c/wrtfmt.c\
	f2c/d_nint.c	  f2c/hl_le.c	    f2c/r_cos.c	      f2c/wsfe.c\
	f2c/d_prod.c	  f2c/hl_lt.c	    f2c/r_cosh.c      f2c/wsle.c\
	f2c/d_sign.c	  f2c/i77vers.c	    f2c/r_dim.c	      f2c/wsne.c\
	f2c/d_sin.c	  f2c/i_abs.c	    f2c/r_exp.c	      f2c/xwsne.c\
	f2c/d_sinh.c	  f2c/i_dim.c	    f2c/r_imag.c      f2c/z_abs.c\
	f2c/d_sqrt.c	  f2c/i_dnnt.c	    f2c/r_int.c	      f2c/z_cos.c\
	f2c/d_tan.c	  f2c/i_indx.c	    f2c/r_lg10.c      f2c/z_div.c\
	f2c/d_tanh.c	  f2c/i_len.c	    f2c/r_log.c	      f2c/z_exp.c\
	f2c/derf_.c	  f2c/i_mod.c	    f2c/r_mod.c	      f2c/z_log.c\
	f2c/derfc_.c	  f2c/i_nint.c	    f2c/r_nint.c      f2c/z_sin.c\
	f2c/dfe.c	  f2c/i_sign.c	    f2c/r_sign.c      f2c/z_sqrt.c

# We also have to pack f2c/arithchk.c in the distribution in case the
# user wants to compile it using the internal f2c. f2c/arithchk.c is
# not linked into libf2c.la (hence we cannot add it to libf2c_la_SOURCES)
# but is needed to build f2c/arith.h
EXTRA_DIST      = f2c/arithchk.c

if INTERNAL_F2C
  libf2c_la_SOURCES = f2c.h f2c/fio.h f2c/fmt.h f2c/sysdep1.h\
					  f2c/sysdep1.h0 f2c/lio.h f2c/fp.h f2c/signal1.h\
					  f2c/signal1.h0 $(F2C)
  libf2c_la_CFLAGS = -DSkip_f2c_Undefs -I. -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
  f2c/arith.h: f2c/arithchk.c
	$(CC) $(CFLAGS) -DNO_FPINIT f2c/arithchk.c -lm -o f2c/arith || \
	$(CC) -DNO_LONG_LONG $(CFLAGS) -DNO_FPINIT f2c/arithchk.c \
		$(WARNING_CFLAGS) -lm -o f2c/arith
	f2c/arith > f2c/arith.h
  f2c/sysdep1.h: f2c/sysdep1.h0
	cp f2c/sysdep1.h0 f2c/sysdep1.h
  f2c/signal1.h: f2c/signal1.h0
	cp f2c/signal1.h0 f2c/signal1.h
  F2C_LIB = libf2c.la

  BUILT_SOURCES += f2c/arith.h
endif

if INTERNAL_BLAS
  libblas_la_SOURCES = f2c.h $(BLAS)
  libblas_la_CFLAGS = -I. -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
  BLAS_LIB = libblas.la
endif

if INTERNAL_LAPACK
  liblapack_la_SOURCES = f2c.h $(LAPACK)
  liblapack_la_CFLAGS = -I. -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
  libdlamch_la_SOURCES = lapack/dlamch.c
  libdlamch_la_CFLAGS = $(FLOATSTORE) -I. -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
  LAPACK_LIB = liblapack.la libdlamch.la
endif

if INTERNAL_ARPACK
  libarpack_la_SOURCES = f2c.h $(ARPACK)
  libarpack_la_CFLAGS = -I. -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
  ARPACK_LIB = libarpack.la
endif

include ../optional/glpk/glpk.inc

if INTERNAL_GLPK
  libglpk_la_SOURCES = 	$(GLPK)
  libglpk_la_CFLAGS = -I$(top_srcdir)/optional/glpk
  libglpk_la_CPPFLAGS = -I$(top_srcdir)/include
  GLPK_LIB = libglpk.la
endif

libplfit_la_SOURCES = $(PLFIT)
PLFIT_LIB = libplfit.la

noinst_LTLIBRARIES = $(F2C_LIB) $(BLAS_LIB) $(LAPACK_LIB) $(ARPACK_LIB) \
					 $(GLPK_LIB) $(PLFIT_LIB)

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
		infomap_Greedy.h infomap_Node.h infomap_Greedy.h infomap_FlowGraph.h \
		igraph_math.h \
		drl_layout.h drl_parse.h drl_graph.h \
		drl_graph_3d.h		drl_layout_3d.h \
		drl_Node.h	drl_Node_3d.h \
		DensityGrid.h	DensityGrid_3d.h \
		igraph_flow_internal.h \
		vector.pmt matrix.pmt stack.pmt dqueue.pmt heap.pmt array.pmt \
		igraph_types_internal.h \
		foreign-dl-header.h bignum.h bigint.h \
		gengraph_box_list.h gengraph_definitions.h \
		gengraph_degree_sequence.h gengraph_graph_molloy_hash.h \
		gengraph_graph_molloy_optimized.h \
		gengraph_hash.h gengraph_header.h gengraph_powerlaw.h \
		gengraph_qsort.h gengraph_random.h gengraph_vertex_cover.h \
		igraph_blas_internal.h igraph_arpack_internal.h \
		igraph_lapack_internal.h igraph_glpk_support.h \
		igraph_marked_queue.h igraph_estack.h \
		hrg_dendro.h hrg_graph.h hrg_rbtree.h hrg_splittree_eq.h \
		hrg_graph_simp.h foreign-gml-header.h \
		foreign-ncol-header.h foreign-lgl-header.h \
		foreign-pajek-header.h igraph_interrupt_internal.h \
		scg_headers.h igraph_hacks_internal.h

HEADERS_PUBLIC =../include/igraph.h 		../include/igraph_memory.h    \
		../include/igraph_random.h 	../include/igraph_types.h     \
		../include/igraph_vector.h 	../include/igraph_matrix.h    \
		../include/igraph_array.h 	../include/igraph_dqueue.h    \
		../include/igraph_stack.h 	../include/igraph_heap.h      \
		../include/igraph_arpack.h 	\
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
		../include/igraph_interface.h	../include/igraph_constants.h \
		../include/igraph_vector_pmt.h	../include/igraph_matrix_pmt.h\
		../include/igraph_array_pmt.h	../include/igraph_dqueue_pmt.h\
		../include/igraph_stack_pmt.h	../include/igraph_heap_pmt.h  \
		../include/igraph_vector_ptr.h	../include/igraph_spmatrix.h  \
		../include/igraph_strvector.h	../include/igraph_psumtree.h  \
		../include/igraph_sparsemat.h   ../include/igraph_mixing.h    \
		../include/igraph_version.h     ../include/igraph_blas.h      \
		../include/igraph_separators.h  ../include/igraph_cohesive_blocks.h \
		../include/igraph_lapack.h	../include/igraph_complex.h   \
		../include/igraph_eigen.h	../include/igraph_statusbar.h \
		../include/igraph_hrg.h         ../include/igraph_microscopic_update.h \
		../include/igraph_interrupt.h   ../include/igraph_threading.h \
		../include/igraph_scg.h		../include/igraph_qsort.h \
		../include/igraph_matching.h

SOURCES = 		     basic_query.c games.c cocitation.c iterators.c \
			     structural_properties.c components.c layout.c \
			     structure_generators.c conversion.c \
			     type_indexededgelist.c spanning_trees.c \
			     igraph_error.c interrupt.c other.c foreign.c random.c \
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
			     infomap.cc infomap_Greedy.cc infomap_Node.cc infomap_FlowGraph.cc \
			     spmatrix.c community.c fast_community.c \
			     evolver_cit.c gml_tree.c \
			     bliss_orbit.cc bliss_eqrefhash.cc \
			     bliss_partition.cc bliss_graph.cc bliss_timer.cc \
			     bliss_heap.cc bliss_utils.cc bliss.cc \
			     cattributes.c revolver_ml_cit.c zeroin.c bfgs.c math.c \
			     forestfire.c microscopic_update.c \
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
			     igraph_fixed_vectorlist.c separators.c \
			     igraph_marked_queue.c igraph_estack.c st-cuts.c \
			     cohesive_blocks.c statusbar.c \
			     lapack.c complex.c eigen.c feedback_arc_set.c \
			     sugiyama.c glpk_support.c \
				 igraph_hrg_types.cc igraph_hrg.cc \
			     distances.c fortran_intrinsics.c matching.c \
			     scg.c scg_approximate_methods.c scg_exact_scg.c \
			     scg_kmeans.c scg_utils.c scg_optimal_method.c \
			     qsort.c qsort_r.c types.c hacks.c

libigraph_la_SOURCES       = $(SOURCES) $(HEADERS_PRIVATE)
libigraph_la_CFLAGS        = -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
libigraph_la_CXXFLAGS	   = -I$(top_srcdir)/include -I$(top_builddir)/include $(WARNING_CFLAGS)
libigraph_la_LDFLAGS       = -no-undefined
libigraph_la_LIBADD        = -lm $(XML2_LIBS) $(F2C_LIB) $(BLAS_LIB) \
				 $(LAPACK_LIB) $(ARPACK_LIB) $(GLPK_LIB) $(PLFIT_LIB)

if INTERNAL_GLPK
  libigraph_la_CFLAGS   += -I$(top_srcdir)/optional/glpk
  libigraph_la_CXXFLAGS += -I$(top_srcdir)/optional/glpk
endif

includedir = $(prefix)/include/igraph
include_HEADERS = $(HEADERS_PUBLIC)

MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h

echosources:
	$(info $(SOURCES) $(LAPACK) $(ARPACK) $(BLAS) $(libf2c_la_SOURCES) \
	$(libdlamch_la_SOURCES) $(libplfit_la_SOURCES))

echoheaders:
	$(info $(HEADERS_PUBLIC))

echoheadersprivate:
	$(info $(HEADERS_PRIVATE))

.PHONY: sources
