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
	lapack/dlarfg.c  lapack/dlasrt.c
  noinst_LTLIBRARIES  = libdlamch.la
  libdlamch_la_SOURCES = lapack/dlamch.c
  libdlamch_la_CFLAGS = -ffloat-store
  LAPACK_DLAMCH = libdlamch.la
else 
  LAPACK =
  LAPACK_DLAMCH = 
endif

if INTERNAL_ARPACK
  ARPACK = arpack/dgetv0.c arpack/dmout.c arpack/dsaitr.c arpack/dsapps.c \
         arpack/dsaup2.c arpack/dsaupd.c arpack/dsconv.c arpack/dseigt.c \
         arpack/dsesrt.c arpack/dseupd.c arpack/dsgets.c arpack/dsortr.c \
         arpack/dstats.c arpack/dstqrb.c arpack/dvout.c arpack/ivout.c \
         arpack/second.c
else
  ARPACK =
endif

if INTERNAL_BLAS
  BLAS = blas/daxpy.c  blas/ddot.c   blas/dger.c   blas/dscal.c  blas/lsame.c \
	blas/dcopy.c  blas/dgemv.c  blas/dnrm2.c  blas/dswap.c
else
  BLAS =
endif

F2C = f2c/d_sign.c f2c/etime_.c  f2c/pow_dd.c f2c/pow_di.c f2c/s_cmp.c \
	f2c/s_copy.c

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
			     arpack.c

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
		  dqueue.pmt dqueue.h heap.pmt heap.h arpack.h

MAINTAINERCLEANFILES = Makefile.in foreign-lgl-parser.h foreign-ncol-parser.h
