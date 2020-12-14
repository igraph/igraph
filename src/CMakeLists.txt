
# Traverse subdirectories
add_subdirectory(bliss)
add_subdirectory(cliquer)
add_subdirectory(cs)
add_subdirectory(f2c)
add_subdirectory(lapack)
add_subdirectory(plfit)
add_subdirectory(prpack)
add_subdirectory(mini-gmp)

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
  core/dqueue.c
  core/error.c
  core/heap.c
  core/indheap.c
  core/interrupt.c
  core/matrix.c
  core/memory.c
  core/printing.c
  core/set.c
  core/stack.c
  core/vector_ptr.c
  core/vector.c

  io/graphml.c
  io/gml-tree.c
  io/io_other.c
  ${PARSER_SOURCES}

  layout/davidson_harel.c
  layout/drl/DensityGrid.cpp
  layout/drl/DensityGrid_3d.cpp
  layout/drl/drl_graph.cpp
  layout/drl/drl_graph_3d.cpp
  layout/drl/drl_layout.cpp
  layout/drl/drl_layout_3d.cpp
  layout/fruchterman_reingold.c
  layout/gem.c
  layout/kamada_kawai.c
  layout/layout_other.c
  layout/sugiyama.c

  math/bfgs.c
  math/bigint.c
  math/bignum.c
  math/complex.c
  math/utils.c

  adjlist.c
  arpack.c
  atlas.c
  attributes.c
  basic_query.c
  bipartite.c
  blas.c
  bliss.cc
  cattributes.c
  centrality.c
  cliques.c
  clustertool.cpp
  cocitation.c
  cohesive_blocks.c
  coloring.c
  community_leiden.c
  community.c
  components.c
  conversion.c
  cores.c
  decomposition.c
  degree_sequence.cpp
  distances.c
  dotproduct.c
  eigen.c
  embedding.c
  fast_community.c
  feedback_arc_set.c
  flow.c
  forestfire.c
  games.c
  gengraph_box_list.cpp
  gengraph_degree_sequence.cpp
  gengraph_graph_molloy_hash.cpp
  gengraph_graph_molloy_optimized.cpp
  gengraph_mr-connected.cpp
  gengraph_powerlaw.cpp
  glet.c
  glpk_support.c
  graphicality.c
  hacks.c
  igraph_buckets.c
  igraph_cliquer.c
  igraph_estack.c
  igraph_eulerian.c
  igraph_fixed_vectorlist.c
  igraph_grid.c
  igraph_hashtable.c
  igraph_hrg_types.cc
  igraph_hrg.cc
  igraph_marked_queue.c
  igraph_psumtree.c
  igraph_strvector.c
  igraph_trie.c
  infomap_FlowGraph.cc
  infomap_Greedy.cc
  infomap_Node.cc
  infomap.cc
  iterators.c
  lad.c
  lapack.c
  lsap.c
  matching.c
  maximal_cliques.c
  microscopic_update.c
  mixing.c
  motifs.c
  NetDataTypes.cpp
  NetRoutines.cpp
  operators.c
  optimal_modularity.c
  other.c
  paths.c
  pottsmodel_2.cpp
  progress.c
  prpack.cpp
  qsort_r.c
  qsort.c
  random_walk.c
  random.c
  sbm.c
  scan.c
  scg_approximate_methods.c
  scg_exact_scg.c
  scg_kmeans.c
  scg_optimal_method.c
  scg_utils.c
  scg.c
  separators.c
  shortest_paths.c
  sir.c
  spanning_trees.c
  sparsemat.c
  spectral_properties.c
  spmatrix.c
  st-cuts.c
  statusbar.c
  structural_properties.c
  structure_generators.c
  topology.c
  triangles.c
  type_indexededgelist.c
  version.c
  visitors.c
  walktrap_communities.cpp
  walktrap_graph.cpp
  walktrap_heap.cpp
  walktrap.cpp
  zeroin.c

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
