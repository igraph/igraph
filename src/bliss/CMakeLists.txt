# Declare the files needed to compile bliss
add_library(
  bliss
  OBJECT
  EXCLUDE_FROM_ALL
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
  ${GMP_INCLUDE_DIR}
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET bliss PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# If we are linking igraph to GMP, then Bliss will use GMP for bignums
if(GMP_FOUND)
  target_link_libraries(
    bliss
    PUBLIC
    ${GMP_LIBRARY}
  )
endif()

use_all_warnings(bliss)
