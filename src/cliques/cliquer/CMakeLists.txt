# Declare the files needed to compile Cliquer
add_library(
  cliquer
  OBJECT
  EXCLUDE_FROM_ALL
  cliquer.c
  cliquer_graph.c
  reorder.c
)

target_include_directories(
  cliquer
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/src
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET cliquer PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(cliquer PRIVATE IGRAPH_STATIC)

# TODO(ntamas): make sure that this works for Cliquer
# use_all_warnings(cliquer)
