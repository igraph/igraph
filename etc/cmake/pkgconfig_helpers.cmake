# Helper functions for generating a nicely formatted igraph.pc file from
# igraph.pc.in

include(JoinPaths)
include(CheckCXXSymbolExists)

# Converts the name of a library file (or framework on macOS) into an
# appropriate linker flag (-lsomething or -framework something.framework).
# Returns the input intact if its extension does not look like a shared or
# static library extension.
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

if(MATH_LIBRARY)
  set(PKGCONFIG_LIBS_PRIVATE "-lm")
else()
  set(PKGCONFIG_LIBS_PRIVATE "")
endif()
set(PKGCONFIG_REQUIRES_PRIVATE "")

if(NOT MSVC)
  check_cxx_symbol_exists(_LIBCPP_VERSION "vector" USING_LIBCXX)
  check_cxx_symbol_exists(__GLIBCXX__ "vector" USING_LIBSTDCXX)
  if(USING_LIBCXX)
    set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lc++")
  elseif(USING_LIBSTDCXX)
    set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lstdc++")
  endif()
endif()

if(IGRAPH_GRAPHML_SUPPORT)
  set(PKGCONFIG_REQUIRES_PRIVATE "${PKGCONFIG_REQUIRES_PRIVATE} libxml-2.0")
endif()
if(NOT IGRAPH_USE_INTERNAL_GMP)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lgmp")
endif()
if(NOT IGRAPH_USE_INTERNAL_BLAS)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lblas")
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
if(NOT IGRAPH_USE_INTERNAL_PLFIT)
  set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -lplfit")
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
