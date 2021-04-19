# Declare the files needed to compile mini-gmp
add_library(
  gmp_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  mini-gmp.c
)

target_include_directories(
  gmp_vendored
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET gmp_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

use_all_warnings(gmp_vendored)

if(MSVC)
  target_compile_options(
    gmp_vendored PRIVATE
    /wd4100  # unreferenced formal parameter
    /wd4127  # conditional expression is constant
    /wd4146  # unary minus operator applied to unsigned type
    /wd4189  # local variable is initialized but not referenced
  )
else()
  target_compile_options(
    gmp_vendored PRIVATE
    $<$<C_COMPILER_ID:GCC,Clang,AppleClang>:-Wno-unused-variable>
  )
endif()

