include(helpers)

include(tls)
include(lto)

option(IGRAPH_GLPK_SUPPORT "Compile igraph with GLPK support" ON)
tristate(IGRAPH_GRAPHML_SUPPORT "Compile igraph with GraphML support" AUTO)
tristate(IGRAPH_OPENMP_SUPPORT "Use OpenMP for parallelization" AUTO)

set(IGRAPH_INTEGER_SIZE AUTO CACHE STRING "Set size of igraph integers")
set_property(CACHE IGRAPH_INTEGER_SIZE PROPERTY STRINGS AUTO 32 64)

if(IGRAPH_INTEGER_SIZE STREQUAL AUTO)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(IGRAPH_INTEGER_SIZE 64)
  else()
    set(IGRAPH_INTEGER_SIZE 32)
  endif()
endif()

option(FLEX_KEEP_LINE_NUMBERS "Keep references to the original line numbers in generated Flex/Bison parser files" OFF)
mark_as_advanced(FLEX_KEEP_LINE_NUMBERS)

option(BUILD_FUZZING "Build fuzz targets and enable fuzzer instrumentation" OFF)
mark_as_advanced(BUILD_FUZZING)
