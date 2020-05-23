# Declare the files needed to compile bluss
add_library(
  bliss
  STATIC
  bliss_heap.cc
  defs.cc
  graph.cc
  orbit.cc
  partition.cc
  uintseqhash.cc
  utils.cc
)

target_include_directories(
  bliss
  PRIVATE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/src
)
