# Declare the files needed to compile bliss
add_library(
  bliss
  OBJECT
  EXCLUDE_FROM_ALL
  defs.cc
  graph.cc
  heap.cc
  orbit.cc
  partition.cc
  uintseqhash.cc
  utils.cc
)

target_include_directories(
  bliss
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/vendor
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/src
  $<$<BOOL:${GMP_INCLUDE_DIR}>:${GMP_INCLUDE_DIR}>
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET bliss PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# If we are linking igraph to GMP, then Bliss will use GMP for bignums
if(NOT GMP_IS_VENDORED)
  target_link_libraries(
    bliss
    PUBLIC
    ${GMP_LIBRARY}
  )
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(bliss PRIVATE IGRAPH_STATIC)

use_all_warnings(bliss)

if (MSVC)
  target_compile_options(bliss PRIVATE /wd4100) # disable unreferenced parameter warning
else()
  target_compile_options(
    bliss PRIVATE
    $<$<C_COMPILER_ID:GCC,Clang,AppleClang>:-Wno-unused-variable>
  )
endif()
