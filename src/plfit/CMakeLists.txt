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

set_property(TARGET plfit PROPERTY POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

use_all_warnings(plfit)
