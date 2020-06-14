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

use_all_warnings(plfit)
