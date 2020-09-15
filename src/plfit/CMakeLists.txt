# Declare the files needed to compile plfit
add_library(
  plfit
  STATIC
  error.c
  gss.c
  hzeta.c
  kolmogorov.c
  lbfgs.c
  mt.c
  options.c
  platform.c
  plfit.c
  rbinom.c
  sampling.c
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET plfit PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

use_all_warnings(plfit)
