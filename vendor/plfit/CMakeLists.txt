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

target_include_directories(
  plfit
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET plfit PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(plfit PRIVATE IGRAPH_STATIC)

use_all_warnings(plfit)

if (MSVC)
  target_compile_options(plfit PRIVATE /wd4100) # disable unreferenced parameter warning
endif()

