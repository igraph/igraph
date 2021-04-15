# Helper functions for generating a nicely formatted igraph.pc file from
# igraph.pc.in

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
    set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} ${OpenMP_${CURRENT_LIB}_LIBRARY}")
  endforeach()
endif()

include(JoinPaths)

function(generate_pkgconfig_file SOURCE TARGET)
  join_paths(PKGCONFIG_LIBDIR "\${exec_prefix}" "${CMAKE_INSTALL_LIBDIR}")
  join_paths(PKGCONFIG_INCLUDEDIR "\${prefix}" "${CMAKE_INSTALL_INCLUDEDIR}")
  configure_file(${SOURCE} ${TARGET})
endfunction()
