# Declare the files needed to compile our vendored copy of the PCG random
# number generator
add_library(
  pcg
  OBJECT
  EXCLUDE_FROM_ALL
  pcg-advance-64.c
  pcg-advance-128.c
  pcg-output-32.c
  pcg-output-64.c
  pcg-output-128.c
  pcg-rngs-64.c
  pcg-rngs-128.c
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET pcg PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

use_all_warnings(pcg)

