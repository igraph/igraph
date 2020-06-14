# Declare the files needed to compile the PRPACK-related stuff
add_library(
  prpack
  STATIC
  prpack_base_graph.cpp
  prpack_igraph_graph.cpp
  prpack_preprocessed_ge_graph.cpp
  prpack_preprocessed_gs_graph.cpp
  prpack_preprocessed_scc_graph.cpp
  prpack_preprocessed_schur_graph.cpp
  prpack_result.cpp
  prpack_solver.cpp
  prpack_utils.cpp
)

target_compile_definitions(
  prpack
  PUBLIC
  PRPACK_IGRAPH_SUPPORT=1
)

target_include_directories(
  prpack
  PRIVATE
  ${CMAKE_SOURCE_DIR}/include
)

# Turn on all warnings for GCC, clang and MSVC
use_all_warnings(prpack)
