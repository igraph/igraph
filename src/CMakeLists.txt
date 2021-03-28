
# Traverse subdirectories
add_subdirectory(centrality/prpack)
add_subdirectory(cliques/cliquer)
add_subdirectory(isomorphism/bliss)

# Generate lexers and parsers
set(PARSER_SOURCES)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/io/parsers)
foreach(FORMAT dl gml lgl ncol pajek)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/io/parsers/${FORMAT}-parser.c)
    list(APPEND PARSER_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/io/parsers/${FORMAT}-lexer.c
      ${CMAKE_CURRENT_SOURCE_DIR}/io/parsers/${FORMAT}-parser.c
    )
  else()
    bison_target(
      ${FORMAT}_parser io/${FORMAT}-parser.y ${CMAKE_CURRENT_BINARY_DIR}/io/parsers/${FORMAT}-parser.c
    )
    flex_target(
      ${FORMAT}_lexer io/${FORMAT}-lexer.l ${CMAKE_CURRENT_BINARY_DIR}/io/parsers/${FORMAT}-lexer.c
    )
    add_flex_bison_dependency(${FORMAT}_lexer ${FORMAT}_parser)
    list(APPEND PARSER_SOURCES ${BISON_${FORMAT}_parser_OUTPUTS} ${FLEX_${FORMAT}_lexer_OUTPUTS})
  endif()
endforeach()
add_custom_target(parsersources SOURCES ${PARSER_SOURCES})

# Declare the files needed to compile the igraph library
add_library(
  igraph
  core/array.c
  core/buckets.c
  core/cutheap.c
  core/dqueue.c
  core/error.c
  core/estack.c
  core/fixed_vectorlist.c
  core/grid.c
  core/hashtable.c
  core/heap.c
  core/indheap.c
  core/interruption.c
  core/marked_queue.c
  core/matrix.c
  core/memory.c
  core/printing.c
  core/progress.c
  core/psumtree.c
  core/set.c
  core/sparsemat.c
  core/spmatrix.c
  core/stack.c
  core/statusbar.c
  core/strvector.c
  core/trie.c
  core/vector_ptr.c
  core/vector.c

  math/bfgs.c
  math/complex.c
  math/utils.c

  linalg/arpack.c
  linalg/blas.c
  linalg/eigen.c
  linalg/lapack.c

  random/random.c

  graph/adjlist.c
  graph/attributes.c
  graph/basic_query.c
  graph/cattributes.c
  graph/iterators.c
  graph/type_indexededgelist.c
  graph/visitors.c

  constructors/adjacency.c
  constructors/atlas.c
  constructors/basic.c
  constructors/de_bruijn.c
  constructors/famous.c
  constructors/full.c
  constructors/kautz.c
  constructors/lcf.c
  constructors/linegraph.c
  constructors/prufer.c
  constructors/regular.c

  games/barabasi.c
  games/callaway_traits.c
  games/citations.c
  games/correlated.c
  games/degree_sequence_vl/gengraph_box_list.cpp
  games/degree_sequence_vl/gengraph_degree_sequence.cpp
  games/degree_sequence_vl/gengraph_graph_molloy_hash.cpp
  games/degree_sequence_vl/gengraph_graph_molloy_optimized.cpp
  games/degree_sequence_vl/gengraph_mr-connected.cpp
  games/degree_sequence_vl/gengraph_powerlaw.cpp
  games/degree_sequence.c
  games/dotproduct.c
  games/erdos_renyi.c
  games/establishment.c
  games/forestfire.c
  games/grg.c
  games/growing_random.c
  games/islands.c
  games/k_regular.c
  games/preference.c
  games/recent_degree.c
  games/sbm.c
  games/static_fitness.c
  games/tree.c
  games/watts_strogatz.c

  centrality/betweenness.c
  centrality/centrality_other.c
  centrality/centralization.c
  centrality/closeness.c
  centrality/coreness.c
  centrality/prpack.cpp

  cliques/cliquer.c
  cliques/cliques.c
  cliques/maximal_cliques.c
  cliques/glet.c

  community/edge_betweenness.c
  community/fast_modularity.c
  community/fluid.c
  community/infomap/infomap_FlowGraph.cc
  community/infomap/infomap_Greedy.cc
  community/infomap/infomap_Node.cc
  community/infomap/infomap.cc
  community/label_propagation.c
  community/leading_eigenvector.c
  community/leiden.c
  community/louvain.c
  community/misc.c
  community/modularity.c
  community/optimal_modularity.c
  community/spinglass/clustertool.cpp
  community/spinglass/NetDataTypes.cpp
  community/spinglass/NetRoutines.cpp
  community/spinglass/pottsmodel_2.cpp
  community/walktrap/walktrap_communities.cpp
  community/walktrap/walktrap_graph.cpp
  community/walktrap/walktrap_heap.cpp
  community/walktrap/walktrap.cpp

  connectivity/cohesive_blocks.c
  connectivity/components.c
  connectivity/separators.c

  flow/flow.c
  flow/st-cuts.c

  hrg/hrg_types.cc
  hrg/hrg.cc

  io/dimacs.c
  io/dl.c
  io/dot.c
  io/edgelist.c
  io/graphml.c
  io/gml-tree.c
  io/gml.c
  io/graphdb.c
  io/leda.c
  io/lgl.c
  io/ncol.c
  io/pajek.c
  ${PARSER_SOURCES}

  layout/bipartite.c
  layout/circular.c
  layout/davidson_harel.c
  layout/drl/DensityGrid.cpp
  layout/drl/DensityGrid_3d.cpp
  layout/drl/drl_graph.cpp
  layout/drl/drl_graph_3d.cpp
  layout/drl/drl_layout.cpp
  layout/drl/drl_layout_3d.cpp
  layout/fruchterman_reingold.c
  layout/gem.c
  layout/graphopt.c
  layout/grid.c
  layout/kamada_kawai.c
  layout/lgl.c
  layout/mds.c
  layout/merge_dla.c
  layout/merge_grid.c
  layout/random.c
  layout/reingold_tilford.c
  layout/sugiyama.c

  operators/add_edge.c
  operators/complementer.c
  operators/compose.c
  operators/connect_neighborhood.c
  operators/contract.c
  operators/difference.c
  operators/disjoint_union.c
  operators/intersection.c
  operators/misc_internal.c
  operators/permute.c
  operators/rewire.c
  operators/rewire_edges.c
  operators/simplify.c
  operators/subgraph.c
  operators/union.c

  paths/all_shortest_paths.c
  paths/bellman_ford.c
  paths/dijkstra.c
  paths/distances.c
  paths/eulerian.c
  paths/histogram.c
  paths/johnson.c
  paths/random_walk.c
  paths/shortest_paths.c
  paths/simple_paths.c
  paths/unweighted.c

  properties/basic.c
  properties/constraint.c
  properties/convergence_degree.c
  properties/dag.c
  properties/degrees.c
  properties/girth.c
  properties/loops.c
  properties/multiplicity.c
  properties/neighborhood.c
  properties/spectral.c
  properties/trees.c
  properties/triangles.c

  scg/scg_approximate_methods.c
  scg/scg_exact_scg.c
  scg/scg_kmeans.c
  scg/scg_optimal_method.c
  scg/scg_utils.c
  scg/scg.c

  isomorphism/bliss.cc
  isomorphism/isoclasses.c
  isomorphism/lad.c
  isomorphism/misc.c
  isomorphism/queries.c
  isomorphism/vf2.c

  misc/bipartite.c
  misc/chordality.c
  misc/cocitation.c
  misc/coloring.c
  misc/conversion.c
  misc/degree_sequence.cpp
  misc/embedding.c
  misc/feedback_arc_set.c
  misc/graphicality.c
  misc/matching.c
  misc/microscopic_update.c
  misc/mixing.c
  misc/motifs.c
  misc/other.c
  misc/scan.c
  misc/sir.c
  misc/spanning_trees.c

  internal/glpk_support.c
  internal/hacks.c
  internal/lsap.c
  internal/qsort_r.c
  internal/qsort.c
  internal/zeroin.c

  version.c

  # Vendored library sources. Yes, this is horrible.
  $<IF:$<OR:$<BOOL:${ARPACK_IS_VENDORED}>,$<BOOL:${BLAS_IS_VENDORED}>,$<BOOL:${LAPACK_IS_VENDORED}>>,$<TARGET_OBJECTS:f2c_vendored>,>
  $<IF:$<BOOL:${ARPACK_IS_VENDORED}>,$<TARGET_OBJECTS:arpack_vendored>,>
  $<IF:$<BOOL:${BLAS_IS_VENDORED}>,$<TARGET_OBJECTS:blas_vendored>,>
  $<IF:$<BOOL:${CXSPARSE_IS_VENDORED}>,$<TARGET_OBJECTS:cxsparse_vendored>,>
  $<IF:$<BOOL:${GLPK_IS_VENDORED}>,$<TARGET_OBJECTS:glpk_vendored>,>
  $<IF:$<BOOL:${GMP_IS_VENDORED}>,$<TARGET_OBJECTS:gmp_vendored>,>
  $<IF:$<BOOL:${LAPACK_IS_VENDORED}>,$<TARGET_OBJECTS:lapack_vendored>,>
)

# Set soname for the library
set_target_properties(igraph PROPERTIES SOVERSION 0)

# Add extra compiler definitions if needed
target_compile_definitions(
  igraph
  PRIVATE
  IGRAPH_VERIFY_FINALLY_STACK=$<IF:$<BOOL:${IGRAPH_VERIFY_FINALLY_STACK}>,1,0>
)

# Make sure that a macro named IGRAPH_FILE_BASENAME is provided in every
# compiler call so we can use these in debug messages without revealing the
# full path of the file on the machine where it was compiled
define_file_basename_for_sources(igraph)

# Add include path. Includes are in ../include but they get installed to
# <prefix>/include/igraph, hence the two options. We also have some private
# includes that are generated at compile time but are not part of the public
# interface.
target_include_directories(
  igraph
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include/igraph>
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/vendor

  # Vendored library include paths
  $<$<BOOL:${CXSPARSE_IS_VENDORED}>:$<TARGET_PROPERTY:cxsparse_vendored,INCLUDE_DIRECTORIES>>
  $<$<BOOL:${GLPK_IS_VENDORED}>:$<TARGET_PROPERTY:glpk_vendored,INCLUDE_DIRECTORIES>>

  # Include paths for dependencies
  $<$<BOOL:${CXSPARSE_INCLUDE_DIRS}>:${CXSPARSE_INCLUDE_DIRS}>
  $<$<BOOL:${GLPK_INCLUDE_DIR}>:${GLPK_INCLUDE_DIR}>
  $<$<BOOL:${GMP_INCLUDE_DIR}>:${GMP_INCLUDE_DIR}>
  $<$<BOOL:${LIBXML2_INCLUDE_DIRS}>:${LIBXML2_INCLUDE_DIRS}>
)

if(MATH_LIBRARY)
  target_link_libraries(igraph PUBLIC ${MATH_LIBRARY})
endif()

if(ARPACK_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${ARPACK_LIBRARIES})
endif()

if(BLAS_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${BLAS_LIBRARIES})
endif()

if(CXSPARSE_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${CXSPARSE_LIBRARIES})
endif()

if(GLPK_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${GLPK_LIBRARIES})
endif()

if(GMP_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${GMP_LIBRARIES})
endif()

if(LAPACK_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${LAPACK_LIBRARIES})
endif()

if(LIBXML2_LIBRARIES)
  target_link_libraries(igraph PUBLIC ${LIBXML2_LIBRARIES})
endif()

# Link igraph statically to some of the libraries from the subdirectories
target_link_libraries(
  igraph
  PRIVATE
  bliss cliquer plfit prpack
)

# Disable complex number support for CXSparse because:
#   - It is necessary to compile with MSVC
#   - igraph does not need complex number support from CXSparse on any platform
# This is needed here (in addition to the cxsparse_vendored target) because
# igraph may be compiled with an external CXSparse.
target_compile_definitions(igraph PRIVATE NCOMPLEX)

if (NOT BUILD_SHARED_LIBS)
  target_compile_definitions(igraph PRIVATE IGRAPH_STATIC)
else()
  target_compile_definitions(igraph PRIVATE igraph_EXPORTS)
endif()

if(MSVC)
  # Add MSVC-specific include path for some headers that are missing on Windows
  target_include_directories(igraph PRIVATE ${PROJECT_SOURCE_DIR}/msvc/include)
endif()

# Turn on all warnings for GCC, clang and MSVC
use_all_warnings(igraph)

# GNUInstallDirs be included before generating the pkgconfig file, as it defines
# CMAKE_INSTALL_LIBDIR and CMAKE_INSTALL_INCLUDEDIR variables.
include(GNUInstallDirs)

# Generate pkgconfig file.
# The library names being used here are Linux-specific, but pkgconfig files
# are omstly used on Linux anyway. Nevertheless, we take care not to include
# -lm on Windows because the Python interface of igraph uses the pkg-config
# file to decide what to link to, and we don't have a separate math library
# on Windows.
if(WIN32)
  set(PKGCONFIG_LIBS_PRIVATE "")
else()
  set(PKGCONFIG_LIBS_PRIVATE "-lm")
endif()
if(APPLE)
  # All recent macOS distributions use libc++
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lc++")
elseif(NOT MSVC)
  # Most Linux distributions and MSYS use libstdc++
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lstdc++")
endif()

if(IGRAPH_GRAPHML_SUPPORT)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lxml2 -lz")
endif()
if(NOT IGRAPH_USE_INTERNAL_GMP)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lgmp")
endif()
if(NOT IGRAPH_USE_INTERNAL_BLAS)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lblas")
endif()
if(NOT IGRAPH_USE_INTERNAL_CXSPARSE)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lcxsparse")
endif()
if(IGRAPH_GLPK_SUPPORT AND NOT IGRAPH_USE_INTERNAL_GLPK)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lglpk")
endif()
if(NOT IGRAPH_USE_INTERNAL_LAPACK)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -llapack")
endif()
if(NOT IGRAPH_USE_INTERNAL_ARPACK)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -larpack")
endif()

include(JoinPaths)
join_paths(PKGCONFIG_LIBDIR "\${exec_prefix}" "${CMAKE_INSTALL_LIBDIR}")
join_paths(PKGCONFIG_INCLUDEDIR "\${prefix}" "${CMAKE_INSTALL_INCLUDEDIR}")
configure_file(${PROJECT_SOURCE_DIR}/igraph.pc.in ${PROJECT_BINARY_DIR}/igraph.pc @ONLY)

include(GenerateExportHeader)
generate_export_header(igraph
  STATIC_DEFINE IGRAPH_STATIC
  EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/igraph_export.h
)

# Define how to install the library
install(
  TARGETS igraph
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/igraph
  FILES_MATCHING PATTERN "*.h"
)
install(
  DIRECTORY ${PROJECT_BINARY_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/igraph
  FILES_MATCHING PATTERN "*.h"
)
install(
  FILES ${PROJECT_BINARY_DIR}/igraph.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
