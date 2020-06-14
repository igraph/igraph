# Include some compiler-related helper functions
include(compilers)

# Traverse subdirectories
add_subdirectory(bliss)
add_subdirectory(cliquer)
add_subdirectory(f2c)
add_subdirectory(lapack)
add_subdirectory(plfit)
add_subdirectory(prpack)

# Generate lexers and parsers
set(PARSER_SOURCES)
foreach(FORMAT dl gml lgl ncol pajek)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/foreign-${FORMAT}-parser.c)
    list(APPEND PARSER_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/foreign-${FORMAT}-lexer.c
      ${CMAKE_CURRENT_SOURCE_DIR}/foreign-${FORMAT}-parser.c
    )
  else()
    bison_target(
      ${FORMAT}_parser foreign-${FORMAT}-parser.y ${CMAKE_CURRENT_BINARY_DIR}/foreign-${FORMAT}-parser.c
    )
    flex_target(
      ${FORMAT}_lexer foreign-${FORMAT}-lexer.l ${CMAKE_CURRENT_BINARY_DIR}/foreign-${FORMAT}-lexer.c
    )
    add_flex_bison_dependency(${FORMAT}_lexer ${FORMAT}_parser)
    list(APPEND PARSER_SOURCES ${BISON_${FORMAT}_parser_OUTPUTS} ${FLEX_${FORMAT}_lexer_OUTPUTS})
  endif()
endforeach()

# Declare the files needed to compile the igraph library
add_library(
  igraph
  adjlist.c
  arpack.c
  array.c
  atlas.c
  attributes.c
  basic_query.c
  bfgs.c
  bigint.c
  bignum.c
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
  complex.c
  components.c
  conversion.c
  cores.c
  decomposition.c
  degree_sequence.cpp
  DensityGrid_3d.cpp
  DensityGrid.cpp
  distances.c
  dotproduct.c
  dqueue.c
  drl_graph.cpp
  drl_graph_3d.cpp
  drl_layout.cpp
  drl_layout_3d.cpp
  eigen.c
  embedding.c
  f2c_dummy.c
  fast_community.c
  feedback_arc_set.c
  flow.c
  foreign-graphml.c
  foreign.c
  forestfire.c
  fortran_intrinsics.c
  games.c
  gengraph_box_list.cpp
  gengraph_degree_sequence.cpp
  gengraph_graph_molloy_hash.cpp
  gengraph_graph_molloy_optimized.cpp
  gengraph_mr-connected.cpp
  gengraph_powerlaw.cpp
  glet.c
  glpk_support.c
  gml_tree.c
  hacks.c
  heap.c
  igraph_buckets.c
  igraph_cliquer.c
  igraph_error.c
  igraph_estack.c
  igraph_fixed_vectorlist.c
  igraph_grid.c
  igraph_hashtable.c
  igraph_heap.c
  igraph_hrg_types.cc
  igraph_hrg.cc
  igraph_marked_queue.c
  igraph_psumtree.c
  igraph_set.c
  igraph_stack.c
  igraph_strvector.c
  igraph_trie.c
  infomap_FlowGraph.cc
  infomap_Greedy.cc
  infomap_Node.cc
  infomap.cc
  interrupt.c
  iterators.c
  lad.c
  lapack.c
  layout_dh.c
  layout_fr.c
  layout_gem.c
  layout_kk.c
  layout.c
  lsap.c
  matching.c
  math.c
  matrix.c
  maximal_cliques.c
  memory.c
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
  sir.c
  spanning_trees.c
  sparsemat.c
  spectral_properties.c
  spmatrix.c
  st-cuts.c
  statusbar.c
  structural_properties.c
  structure_generators.c
  sugiyama.c
  topology.c
  triangles.c
  type_indexededgelist.c
  types.c
  vector_ptr.c
  vector.c
  version.c
  visitors.c
  walktrap_communities.cpp
  walktrap_graph.cpp
  walktrap_heap.cpp
  walktrap.cpp
  zeroin.c

  ${PARSER_SOURCES}
)

# Set soname for the library
set_target_properties(igraph PROPERTIES SOVERSION 0)

# Add extra compiler definitions if needed
target_compile_definitions(
  igraph
  PRIVATE
  IGRAPH_VERIFY_FINALLY_STACK=$<IF:$<BOOL:${IGRAPH_VERIFY_FINALLY_STACK}>,1,0>
)

# Add include path. Includes are in ../include but they get installed to
# <prefix/include/igraph, hence the two options. We also have some private
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
  ${CXSPARSE_INCLUDE_DIRECTORIES}
  ${GLPK_INCLUDE_DIR}
  ${GMP_INCLUDE_DIR}
  ${LIBXML2_INCLUDE_DIRS}
)

# Link igraph statically to some of the libraries from the subdirectories
target_link_libraries(
  igraph
  PUBLIC
  ${ARPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${CXSPARSE_LIBRARIES}
  ${GLPK_LIBRARIES}
  ${GMP_LIBRARY}
  ${LAPACK_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  PRIVATE
  bliss cliquer plfit prpack
)

# Turn on all warnings for GCC, clang and MSVC
use_all_warnings(igraph)

# Generate pkgconfig file
configure_file(${CMAKE_SOURCE_DIR}/igraph.pc.cmake.in ${CMAKE_BINARY_DIR}/igraph.pc @ONLY)

# Define how to install the library
include(GNUInstallDirs)
install(
  TARGETS igraph
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)
install(
  DIRECTORY ${CMAKE_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
)
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
)
install(
  FILES ${CMAKE_BINARY_DIR}/igraph.pc
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig
)
