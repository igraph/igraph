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

# Check if GCC-style enum value deprecation is supported

include(CheckCSourceCompiles)

check_c_source_compiles(
  "enum { A __attribute__ ((deprecated)) = 0 }; int main() { return 0; }"
  ENUMVAL_DEPRECATION_SUPPORTED
)

if(ENUMVAL_DEPRECATION_SUPPORTED)
  set(IGRAPH_DEPRECATED_ENUMVAL "__attribute__ ((deprecated))")
endif()
