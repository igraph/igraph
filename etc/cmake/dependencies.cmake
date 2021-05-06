include(helpers)

include(CheckSymbolExists)

# The threading library is not needed for igraph itself, but might be needed
# for tests
include(FindThreads)

macro(find_dependencies)
  # Declare the list of dependencies that _may_ be vendored and those that may not
  set(VENDORABLE_DEPENDENCIES BLAS CXSparse GLPK LAPACK ARPACK GMP)
  set(NONVENDORABLE_DEPENDENCIES GLPK OpenMP)

  # Declare configuration options for dependencies
  tristate(IGRAPH_USE_INTERNAL_GMP "Compile igraph with internal Mini-GMP" AUTO)
  tristate(IGRAPH_USE_INTERNAL_ARPACK "Compile igraph with internal ARPACK" AUTO)
  tristate(IGRAPH_USE_INTERNAL_BLAS "Compile igraph with internal BLAS" AUTO)
  tristate(IGRAPH_USE_INTERNAL_CXSPARSE "Compile igraph with internal CXSparse" AUTO)
  tristate(IGRAPH_USE_INTERNAL_GLPK "Compile igraph with internal GLPK" AUTO)
  tristate(IGRAPH_USE_INTERNAL_LAPACK "Compile igraph with internal LAPACK" AUTO)

  # Declare dependencies
  set(REQUIRED_DEPENDENCIES "")
  set(OPTIONAL_DEPENDENCIES FLEX BISON OpenMP)
  set(VENDORED_DEPENDENCIES "")

  # Extend dependencies depending on whether we will be using the vendored
  # copies or not
  foreach(DEPENDENCY ${VENDORABLE_DEPENDENCIES})
    string(TOUPPER "${DEPENDENCY}" LIBNAME_UPPER)

    if(IGRAPH_USE_INTERNAL_${LIBNAME_UPPER} STREQUAL "AUTO")
      find_package(${DEPENDENCY})
      if(${LIBNAME_UPPER}_FOUND)
        set(IGRAPH_USE_INTERNAL_${LIBNAME_UPPER} OFF)
      else()
        set(IGRAPH_USE_INTERNAL_${LIBNAME_UPPER} ON)
      endif()
    endif()

    if(IGRAPH_USE_INTERNAL_${LIBNAME_UPPER})
      list(APPEND VENDORED_DEPENDENCIES ${DEPENDENCY})
    else()
      list(APPEND REQUIRED_DEPENDENCIES ${DEPENDENCY})
    endif()
  endforeach()

  # For nonvendorable dependencies, figure out whether we should attempt to
  # link to them based on the value of the IGRAPH_..._SUPPORT option
  foreach(DEPENDENCY ${NONVENDORABLE_DEPENDENCIES})
    string(TOUPPER "${DEPENDENCY}" LIBNAME_UPPER)

    if(IGRAPH_${LIBNAME_UPPER}_SUPPORT STREQUAL "AUTO")
      find_package(${DEPENDENCY})
      if(${LIBNAME_UPPER}_FOUND)
        set(IGRAPH_${LIBNAME_UPPER}_SUPPORT ON)
      else()
        set(IGRAPH_${LIBNAME_UPPER}_SUPPORT OFF)
      endif()
    endif()
  endforeach()

  # GraphML support is treated separately because the library name is different
  if(IGRAPH_GRAPHML_SUPPORT STREQUAL "AUTO")
    find_package(LibXml2)
    if(LibXml2_FOUND)
      set(IGRAPH_GRAPHML_SUPPORT ON)
    else()
      set(IGRAPH_GRAPHML_SUPPORT OFF)
    endif()
  endif()

  if(NOT IGRAPH_GLPK_SUPPORT)
    if(IGRAPH_USE_INTERNAL_GLPK)
      list(REMOVE_ITEM VENDORED_DEPENDENCIES GLPK)
    else()
      list(REMOVE_ITEM REQUIRED_DEPENDENCIES GLPK)
    endif()
  endif()

  if(IGRAPH_GRAPHML_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES LibXml2)
  endif()

  # Find dependencies
  foreach(DEPENDENCY ${REQUIRED_DEPENDENCIES} ${OPTIONAL_DEPENDENCIES})
    list(FIND REQUIRED_DEPENDENCIES "${DEPENDENCY}" INDEX)
    set(NEED_THIS_DEPENDENCY NO)

    if(INDEX GREATER_EQUAL 0)
      # This is a required dependency, search for it unconditionally. Do
      # not use REQUIRED; we will report errors in a single batch at the end
      # of the configuration process
      set(NEED_THIS_DEPENDENCY YES)
    else()
      # This is an optional dependency, search for it only if the user did not
      # turn it off explicitly
      string(TOUPPER "${DEPENDENCY}" LIBNAME_UPPER)
      if(NOT DEFINED IGRAPH_${LIBNAME_UPPER}_SUPPORT)
        set(NEED_THIS_DEPENDENCY YES)
      elseif(IGRAPH_${LIBNAME_UPPER}_SUPPORT)
        set(NEED_THIS_DEPENDENCY YES)
      endif()
    endif()

    if(NEED_THIS_DEPENDENCY AND NOT DEFINED ${DEPENDENCY}_FOUND)
      find_package(${DEPENDENCY})
    endif()
  endforeach()

  # Override libraries of vendored dependencies even if they were somehow
  # detected above
  foreach(DEPENDENCY ${VENDORED_DEPENDENCIES})
    string(TOUPPER "${DEPENDENCY}" LIBNAME_UPPER)
    string(TOLOWER "${DEPENDENCY}" LIBNAME_LOWER)
    if(IGRAPH_USE_INTERNAL_${LIBNAME_UPPER})
      set(${LIBNAME_UPPER}_LIBRARIES "")
      set(${LIBNAME_UPPER}_FOUND 1)
      set(${LIBNAME_UPPER}_IS_VENDORED 1)
      set(INTERNAL_${LIBNAME_UPPER} 1)
    endif()
  endforeach()

  # Export some aliases that will be used in config.h
  set(HAVE_GLPK ${GLPK_FOUND})
  set(HAVE_GMP ${GMP_FOUND})
  set(HAVE_LIBXML ${LIBXML2_FOUND})

  # Check whether we need to link to the math library
  if(NOT DEFINED CACHE{NEED_LINKING_AGAINST_LIBM})
    set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
    set(CMAKE_REQUIRED_QUIET ON)
    check_symbol_exists(sinh "math.h" SINH_FUNCTION_EXISTS)
    if(NOT SINH_FUNCTION_EXISTS)
      unset(SINH_FUNCTION_EXISTS CACHE)
      set(CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
      list(APPEND CMAKE_REQUIRED_LIBRARIES m)
      check_symbol_exists(sinh "math.h" SINH_FUNCTION_EXISTS)
      if(SINH_FUNCTION_EXISTS)
        set(NEED_LINKING_AGAINST_LIBM True CACHE BOOL "" FORCE)
      else()
        message(FATAL_ERROR "Failed to figure out how to link to the math library on this platform")
      endif()
      set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})
    endif()
    unset(SINH_FUNCTION_EXISTS CACHE)
	set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})
  endif()
  
  if(NEED_LINKING_AGAINST_LIBM)
    find_library(MATH_LIBRARY m)
  endif()

  mark_as_advanced(MATH_LIBRARY)
  mark_as_advanced(NEED_LINKING_AGAINST_LIBM)
endmacro()
