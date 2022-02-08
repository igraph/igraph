# Declare the files needed to compile our vendored plfit copy
add_library(
  plfit_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  gss.c
  hzeta.c
  kolmogorov.c
  lbfgs.c
  mt.c
  options.c
  platform.c
  plfit.c
  plfit_error.c
  rbinom.c
  sampling.c
)

target_include_directories(
  plfit_vendored
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET plfit_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(plfit_vendored PRIVATE IGRAPH_STATIC)

use_all_warnings(plfit_vendored)

if (MSVC)
  target_compile_options(
    plfit_vendored PRIVATE
    /wd4100
  ) # disable unreferenced parameter warning
endif()

if(IGRAPH_OPENMP_SUPPORT)
  target_link_libraries(plfit_vendored PRIVATE OpenMP::OpenMP_C)
endif()
