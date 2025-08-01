
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
    if (BISON_VERSION VERSION_GREATER_EQUAL 3)
      set(bison_no_deprecated -Wno-deprecated)
    endif()
    if (NOT FLEX_KEEP_LINE_NUMBERS)
      set(bison_hide_line_numbers --no-lines)
      set(flex_hide_line_numbers --noline)
    endif()
    bison_target(
      ${FORMAT}_parser io/${FORMAT}-parser.y ${CMAKE_CURRENT_BINARY_DIR}/io/parsers/${FORMAT}-parser.c
      COMPILE_FLAGS "${bison_hide_line_numbers} ${bison_no_deprecated}"
    )
    flex_target(
      ${FORMAT}_lexer io/${FORMAT}-lexer.l ${CMAKE_CURRENT_BINARY_DIR}/io/parsers/${FORMAT}-lexer.c
      COMPILE_FLAGS "${flex_hide_line_numbers}"
      DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/io/parsers/${FORMAT}-lexer.h
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
  core/bitset.c
  core/bitset_list.c
  core/buckets.c
  core/cutheap.c
  core/dqueue.c
  core/error.c
  core/estack.c
  core/fixed_vectorlist.c
  core/genheap.c
  core/grid.c
  core/heap.c
  core/indheap.c
  core/interruption.c
  core/marked_queue.c
  core/matrix.c
  core/matrix_list.c
  core/memory.c
  core/printing.c
  core/progress.c
  core/psumtree.c
  core/set.c
  core/sparsemat.c
  core/stack.c
  core/statusbar.c
  core/strvector.c
  core/trie.c
  core/vector.c
  core/vector_list.c
  core/vector_ptr.c

  math/complex.c
  math/safe_intop.c
  math/utils.c

  linalg/arpack.c
  linalg/blas.c
  linalg/eigen.c
  linalg/lapack.c

  random/random.c
  random/rng_glibc2.c
  random/rng_mt19937.c
  random/rng_pcg32.c
  random/rng_pcg64.c

  graph/adjlist.c
  graph/attributes.c
  graph/basic_query.c
  graph/caching.c
  graph/cattributes.c
  graph/graph_list.c
  graph/iterators.c
  graph/type_common.c
  graph/type_indexededgelist.c
  graph/visitors.c

  constructors/adjacency.c
  constructors/atlas.c
  constructors/basic_constructors.c
  constructors/circulant.c
  constructors/de_bruijn.c
  constructors/famous.c
  constructors/full.c
  constructors/generalized_petersen.c
  constructors/kautz.c
  constructors/lattices.c
  constructors/lcf.c
  constructors/linegraph.c
  constructors/mycielskian.c
  constructors/prufer.c
  constructors/regular.c
  constructors/trees.c


  games/barabasi.c
  games/callaway_traits.c
  games/chung_lu.c
  games/citations.c
  games/correlated.c
  games/degree_sequence_vl/gengraph_degree_sequence.cpp
  games/degree_sequence_vl/gengraph_graph_molloy_hash.cpp
  games/degree_sequence_vl/gengraph_graph_molloy_optimized.cpp
  games/degree_sequence_vl/gengraph_mr-connected.cpp
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
  centrality/eigenvector.c
  centrality/hub_authority.c
  centrality/pagerank.c
  centrality/truss.cpp
  centrality/prpack.cpp

  cliques/cliquer_wrapper.c
  cliques/cliques.c
  cliques/maximal_cliques.c
  cliques/glet.c

  community/community_misc.c
  community/edge_betweenness.c
  community/fast_modularity.c
  community/fluid.c
  community/infomap/infomap_FlowGraph.cc
  community/infomap/infomap_Greedy.cc
  community/infomap/infomap.cc
  community/label_propagation.c
  community/leading_eigenvector.c
  community/leiden.c
  community/louvain.c
  community/modularity.c
  community/optimal_modularity.c
  community/spinglass/clustertool.cpp
  community/spinglass/NetDataTypes.cpp
  community/spinglass/NetRoutines.cpp
  community/spinglass/pottsmodel_2.cpp
  community/voronoi.c
  community/walktrap/walktrap_communities.cpp
  community/walktrap/walktrap_graph.cpp
  community/walktrap/walktrap_heap.cpp
  community/walktrap/walktrap.cpp

  connectivity/cohesive_blocks.c
  connectivity/components.c
  connectivity/percolation.c
  connectivity/separators.c
  connectivity/reachability.c

  cycles/simple_cycles.c

  flow/flow.c
  flow/flow_conversion.c
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
  io/parse_utils.c
  ${PARSER_SOURCES}

  layout/align.c
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
  layout/kamada_kawai.c
  layout/large_graph.c
  layout/layout_bipartite.c
  layout/layout_grid.c
  layout/layout_random.c
  layout/mds.c
  layout/merge_dla.c
  layout/merge_grid.c
  layout/reingold_tilford.c
  layout/sugiyama.c
  layout/umap.c

  operators/add_edge.c
  operators/complementer.c
  operators/compose.c
  operators/connect_neighborhood.c
  operators/contract.c
  operators/difference.c
  operators/disjoint_union.c
  operators/intersection.c
  operators/join.c
  operators/misc_internal.c
  operators/permute.c
  operators/products.c
  operators/reverse.c
  operators/rewire.c
  operators/rewire_edges.c
  operators/simplify.c
  operators/subgraph.c
  operators/union.c

  paths/all_shortest_paths.c
  paths/astar.c
  paths/bellman_ford.c
  paths/dijkstra.c
  paths/distances.c
  paths/eulerian.c
  paths/floyd_warshall.c
  paths/histogram.c
  paths/johnson.c
  paths/random_walk.c
  paths/shortest_paths.c
  paths/simple_paths.c
  paths/sparsifier.c
  paths/unweighted.c
  paths/voronoi.c
  paths/widest_paths.c

  properties/basic_properties.c
  properties/complete.c
  properties/constraint.c
  properties/convergence_degree.c
  properties/dag.c
  properties/degrees.c
  properties/ecc.c
  properties/girth.c
  properties/loops.c
  properties/multiplicity.c
  properties/neighborhood.c
  properties/perfect.c
  properties/rich_club.c
  properties/spectral.c
  properties/trees.c
  properties/triangles.c

  isomorphism/bliss.cc
  isomorphism/isoclasses.c
  isomorphism/lad.c
  isomorphism/isomorphism_misc.c
  isomorphism/queries.c
  isomorphism/vf2.c

  misc/bipartite.c
  misc/chordality.c
  misc/cocitation.c
  misc/coloring.c
  misc/conversion.c
  misc/cycle_bases.c
  misc/degree_sequence.cpp
  misc/embedding.c
  misc/feedback_arc_set.c
  misc/graphicality.c
  misc/matching.c
  misc/microscopic_update.c
  misc/mixing.c
  misc/motifs.c
  misc/order_cycle.cpp
  misc/other.c
  misc/power_law_fit.c
  misc/scan.c
  misc/sir.c
  misc/spanning_trees.c

  internal/glpk_support.c
  internal/hacks.c
  internal/lsap.c
  internal/qsort_r.c
  internal/qsort.c
  internal/utils.c
  internal/zeroin.c

  version.c

  # Vendored library sources. Yes, this is horrible.
  $<IF:$<OR:$<BOOL:${ARPACK_IS_VENDORED}>,$<BOOL:${BLAS_IS_VENDORED}>,$<BOOL:${LAPACK_IS_VENDORED}>>,$<TARGET_OBJECTS:f2c_vendored>,>
  $<IF:$<BOOL:${ARPACK_IS_VENDORED}>,$<TARGET_OBJECTS:arpack_vendored>,>
  $<IF:$<BOOL:${BLAS_IS_VENDORED}>,$<TARGET_OBJECTS:blas_vendored>,>
  $<IF:$<BOOL:${GLPK_IS_VENDORED}>,$<TARGET_OBJECTS:glpk_vendored>,>
  $<IF:$<BOOL:${GMP_IS_VENDORED}>,$<TARGET_OBJECTS:gmp_vendored>,>
  $<IF:$<BOOL:${LAPACK_IS_VENDORED}>,$<TARGET_OBJECTS:lapack_vendored>,>
  $<IF:$<BOOL:${PLFIT_IS_VENDORED}>,$<TARGET_OBJECTS:plfit_vendored>,>
)

# Required by Xcode new build system
add_dependencies(igraph parsersources)

# Set soname for the library
set_target_properties(igraph PROPERTIES VERSION "3.1.10")
set_target_properties(igraph PROPERTIES SOVERSION 3)

# Add extra compiler definitions if needed
target_compile_definitions(
  igraph
  PRIVATE
  IGRAPH_VERIFY_FINALLY_STACK=$<IF:$<BOOL:${IGRAPH_VERIFY_FINALLY_STACK}>,1,0>
)
# target_compile_options(
#   # -Wconversion could be useful?
#   igraph PRIVATE -Wshorten-64-to-32
# )

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
  "$<$<BOOL:${GLPK_IS_VENDORED}>:$<TARGET_PROPERTY:glpk_vendored,INCLUDE_DIRECTORIES>>"
  "$<$<BOOL:${PLFIT_IS_VENDORED}>:$<TARGET_PROPERTY:plfit_vendored,INCLUDE_DIRECTORIES>>"

  # Include paths for dependencies
  "$<$<BOOL:${GLPK_INCLUDE_DIR}>:${GLPK_INCLUDE_DIR}>"
  "$<$<BOOL:${GMP_INCLUDE_DIR}>:${GMP_INCLUDE_DIR}>"
  "$<$<BOOL:${LIBXML2_INCLUDE_DIRS}>:${LIBXML2_INCLUDE_DIRS}>"
  "$<$<BOOL:${PLFIT_INCLUDE_DIRS}>:${PLFIT_INCLUDE_DIRS}>"
)

if(MATH_LIBRARY)
  target_link_libraries(igraph PUBLIC ${MATH_LIBRARY})
endif()

if(ARPACK_LIBRARIES)
  target_link_libraries(igraph PRIVATE ${ARPACK_LIBRARIES})
endif()

if(BLAS_FOUND)
  target_link_libraries(igraph PRIVATE ${BLAS_LIBRARIES})
endif()

if(GLPK_LIBRARIES)
  target_link_libraries(igraph PRIVATE ${GLPK_LIBRARIES})
endif()

if(GMP_LIBRARIES)
  target_link_libraries(igraph PRIVATE ${GMP_LIBRARIES})
endif()

if(LAPACK_LIBRARIES)
  target_link_libraries(igraph PRIVATE ${LAPACK_LIBRARIES})
endif()

if(LIBXML2_LIBRARIES)
  target_link_libraries(igraph PRIVATE ${LIBXML2_LIBRARIES})
endif()

if(PLFIT_LIBRARIES)
  target_link_libraries(igraph PRIVATE ${PLFIT_LIBRARIES})
endif()

# Link igraph statically to some of the libraries from the subdirectories
target_link_libraries(
  igraph
  PRIVATE
  bliss cliquer cxsparse_vendored pcg prpack
)

if (NOT BUILD_SHARED_LIBS)
  target_compile_definitions(igraph PUBLIC IGRAPH_STATIC)
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

# Generate pkgconfig file
include(pkgconfig_helpers)

include(GenerateExportHeader)
generate_export_header(igraph
  STATIC_DEFINE IGRAPH_STATIC
  EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/igraph_export.h
)

# Provide an igraph-config.cmake file in the installation directory so
# users can find the installed igraph library with FIND_PACKAGE(igraph)
# from their CMakeLists.txt files
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/etc/cmake/igraph-config.cmake.in
  ${PROJECT_BINARY_DIR}/igraph-config.cmake
  # Install destination selected according to https://wiki.debian.org/CMake
  # and by looking at how eigen3 does it in Ubuntu
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/igraph
)
write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/igraph-config-version.cmake
  VERSION ${PACKAGE_VERSION_BASE}
  COMPATIBILITY SameMinorVersion
)

# Define how to install the library
install(
  TARGETS igraph bliss cliquer cxsparse_vendored pcg prpack
  EXPORT igraph_targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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
install(
  FILES ${PROJECT_BINARY_DIR}/igraph-config.cmake
        ${PROJECT_BINARY_DIR}/igraph-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/igraph
)
install(
  EXPORT igraph_targets
  FILE igraph-targets.cmake
  NAMESPACE igraph::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/igraph
)
