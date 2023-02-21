
# Detect if certain attributes are supported by the compiler
# The result will be used to set macros in include/igraph_config.h

# GCC-style enum value deprecation

include(CheckCSourceCompiles)
include(CMakePushCheckState)

# Only check with Clang and GCC as we assume that the -Werror option is supported
# For other compilers, assume that the attribute is unsupported.
if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
  cmake_push_check_state()
  # Require compiling with no warning:
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -Werror")
  check_c_source_compiles(
    "enum { A __attribute__ ((deprecated)) = 0 }; int main(void) { return 0; }"
    COMPILER_HAS_DEPRECATED_ENUMVAL_ATTR
  )
  cmake_pop_check_state()
else()
  set(COMPILER_HAS_DEPRECATED_ENUMVAL_ATTR FALSE)
endif()

if(COMPILER_HAS_DEPRECATED_ENUMVAL_ATTR)
  set(IGRAPH_DEPRECATED_ENUMVAL "__attribute__ ((deprecated))")
endif()
