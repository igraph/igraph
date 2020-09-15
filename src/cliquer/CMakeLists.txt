# Declare the files needed to compile Cliquer
add_library(
  cliquer
  STATIC
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

if (BUILD_SHARED_LIBS)
  set_property(TARGET cliquer PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# TODO(ntamas): make sure that this works for Cliquer
# use_all_warnings(cliquer)
