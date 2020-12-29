
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

  constructors/atlas.c
  constructors/constructors_other.c

  games/degree_sequence_vl/gengraph_box_list.cpp
  games/degree_sequence_vl/gengraph_degree_sequence.cpp
  games/degree_sequence_vl/gengraph_graph_molloy_hash.cpp
  games/degree_sequence_vl/gengraph_graph_molloy_optimized.cpp
  games/degree_sequence_vl/gengraph_mr-connected.cpp
  games/degree_sequence_vl/gengraph_powerlaw.cpp
  games/dotproduct.c
  games/forestfire.c
  games/games_other.c
  games/sbm.c

  centrality/centrality_other.c
  centrality/coreness.c
  centrality/prpack.cpp

  cliques/cliquer.c
  cliques/cliques.c
  cliques/maximal_cliques.c
  cliques/glet.c

  community/community_other.c
  community/fast_modularity.c
  community/leiden.c
  community/optimal_modularity.c
  community/infomap/infomap_FlowGraph.cc
  community/infomap/infomap_Greedy.cc
  community/infomap/infomap_Node.cc
  community/infomap/infomap.cc
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
  layout/layout_other.c
  layout/lgl.c
  layout/mds.c
  layout/merge_dla.c
  layout/random.c
  layout/reingold_tilford.c
  layout/sugiyama.c

  # TODO: each operator should be in a separate file, they are complex enough as they are
  operators/operators.c

  paths/distances.c
  paths/eulerian.c
  paths/random_walk.c
  paths/shortest_paths.c
  paths/simple_paths.c

  properties/properties_other.c
  properties/spectral.c
  properties/triangles.c

  scg/scg_approximate_methods.c
  scg/scg_exact_scg.c
  scg/scg_kmeans.c
  scg/scg_optimal_method.c
  scg/scg_utils.c
  scg/scg.c

  isomorphism/bliss.cc
  isomorphism/lad.c
  isomorphism/isomorphism_other.c

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
  $<TARGET_PROPERTY:prpack,COMPILE_DEFINITIONS>
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
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include/igraph>
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/vendor

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

if (WIN32)
  # Define symbols to enable compilation of CSXSparse on Windows
  target_compile_definitions(igraph PRIVATE NCOMPLEX)
  if (NOT BUILD_SHARED_LIBS)
    # Add a compiler definition required to compile igraph in static mode on Windows
    target_compile_definitions(igraph PRIVATE IGRAPH_STATIC)
  else()
    target_compile_definitions(igraph PRIVATE IGRAPH_EXPORTS)
  endif()
endif()

if(MSVC)
  # Add MSVC-specific include path for some headers that are missing on Windows
  target_include_directories(igraph PRIVATE ${CMAKE_SOURCE_DIR}/msvc/include)
endif()

# Turn on all warnings for GCC, clang and MSVC
use_all_warnings(igraph)

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
configure_file(${CMAKE_SOURCE_DIR}/igraph.pc.in ${CMAKE_BINARY_DIR}/igraph.pc @ONLY)

# Define how to install the library
include(GNUInstallDirs)
install(
  TARGETS igraph
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/igraph
  FILES_MATCHING PATTERN "*.h"
)
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/igraph
  FILES_MATCHING PATTERN "*.h"
)
install(
  FILES ${CMAKE_BINARY_DIR}/igraph.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
