# Declare the files needed to compile plfit
add_library(
  plfit
  OBJECT
  EXCLUDE_FROM_ALL
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

if (MSVC)
  target_compile_options(plfit PRIVATE /wd4100) # disable unreferenced parameter warning
endif()

use_all_warnings(plfit)
