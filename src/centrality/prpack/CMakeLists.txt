# Declare the files needed to compile the PRPACK-related stuff
add_library(
  prpack
  OBJECT
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
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET prpack PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(prpack PRIVATE IGRAPH_STATIC)

# Turn on all warnings for GCC, clang and MSVC
# TODO(ntamas): it does not work yet, gcc errors out with:
#     ignoring #pragma omp parallel [-Werror=unknown-pragmas]
# use_all_warnings(prpack)
