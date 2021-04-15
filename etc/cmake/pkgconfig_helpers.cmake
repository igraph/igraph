# Helper functions for generating a nicely formatted igraph.pc file from
# igraph.pc.in

include(JoinPaths)

# Converts the name of a library file (or framework on macOS) into an
# appropriate linker flag if the library file resides in a standard
# library path. Returns the input intact otherwise.
function(convert_library_file_to_flags output_variable input)
  get_filename_component(input_filename ${input} NAME_WE)
  get_filename_component(input_extension ${input} LAST_EXT)
  if(input_extension STREQUAL ${CMAKE_SHARED_LIBRARY_SUFFIX} OR input_extension STREQUAL ${CMAKE_STATIC_LIBRARY_SUFFIX})
    string(REGEX REPLACE "^${CMAKE_SHARED_LIBRARY_PREFIX}" "" input_stripped ${input_filename})
    set("${output_variable}" "-l${input_stripped}" PARENT_SCOPE)
  elseif(APPLE AND input_extension STREQUAL ".framework")
    set("${output_variable}" "-framework ${input_filename}" PARENT_SCOPE)
  else()
    set("${output_variable}" "${input}" PARENT_SCOPE)
  endif()
endfunction()

# The library names being used here are Linux-specific, but pkgconfig files
# are mostly used on Linux anyway. Nevertheless, we take care not to include
# -lm on Windows because the Python interface of igraph uses the pkg-config
# file to decide what to link to, and we don't have a separate math library
# on Windows.
if(WIN32)
  set(PKGCONFIG_LIBS_PRIVATE "")
else()
  set(PKGCONFIG_LIBS_PRIVATE "-lm")
endif()

if(APPLE)
  # All recent macOS distributions use libc++
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lc++")
elseif(NOT MSVC)
  # Most Linux distributions and MSYS use libstdc++
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lstdc++")
endif()

if(IGRAPH_GRAPHML_SUPPORT)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lxml2 -lz")
endif()
if(NOT IGRAPH_USE_INTERNAL_GMP)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lgmp")
endif()
if(NOT IGRAPH_USE_INTERNAL_BLAS)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lblas")
endif()
if(NOT IGRAPH_USE_INTERNAL_CXSPARSE)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lcxsparse")
endif()
if(IGRAPH_GLPK_SUPPORT AND NOT IGRAPH_USE_INTERNAL_GLPK)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lglpk")
endif()
if(NOT IGRAPH_USE_INTERNAL_LAPACK)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -llapack")
endif()
if(NOT IGRAPH_USE_INTERNAL_ARPACK)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -larpack")
endif()
if(IGRAPH_OPENMP_SUPPORT AND OpenMP_FOUND)
  foreach(CURRENT_LIB ${OpenMP_C_LIB_NAMES})
    convert_library_file_to_flags(CURRENT_LIB "${OpenMP_${CURRENT_LIB}_LIBRARY}")
    set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} ${CURRENT_LIB}")
  endforeach()
endif()

join_paths(PKGCONFIG_LIBDIR "\${exec_prefix}" "${CMAKE_INSTALL_LIBDIR}")
join_paths(PKGCONFIG_INCLUDEDIR "\${prefix}" "${CMAKE_INSTALL_INCLUDEDIR}")
configure_file(
  ${PROJECT_SOURCE_DIR}/igraph.pc.in
  ${PROJECT_BINARY_DIR}/igraph.pc
  @ONLY
)
