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
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/src
)

set_property(TARGET cliquer PROPERTY POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

# TODO(ntamas): make sure that this works for Cliquer
# use_all_warnings(cliquer)
