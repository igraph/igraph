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

use_all_warnings(gmp_vendored)
