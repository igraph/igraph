include(PadString)

macro(find_dependencies)
  # Declare the list of dependencies that _may_ be vendored
  set(VENDORABLE_DEPENDENCIES BLAS LAPACK ARPACK)

  # Declare configuration options for dependencies
  option(IGRAPH_GLPK_SUPPORT "Compile igraph with GLPK support" YES)
  option(IGRAPH_GMP_SUPPORT "Compile igraph with GMP support" YES)
  option(IGRAPH_GRAPHML_SUPPORT "Compile igraph with GraphML support" YES)
  option(IGRAPH_USE_INTERNAL_BLAS "Compile igraph with internal BLAS" NO)
  option(IGRAPH_USE_INTERNAL_LAPACK "Compile igraph with internal LAPACK" NO)
  option(IGRAPH_USE_INTERNAL_ARPACK "Compile igraph with internal ARPACK" NO)

  # Declare dependencies
  set(REQUIRED_DEPENDENCIES CXSparse)
  set(OPTIONAL_DEPENDENCIES FLEX BISON)
  set(VENDORED_DEPENDENCIES "")

  # Extend dependencies depending on whether we will be using the vendored
  # copies or not
  foreach(DEPENDENCY ${VENDORABLE_DEPENDENCIES})
    if(IGRAPH_USE_INTERNAL_${DEPENDENCY})
      list(APPEND VENDORED_DEPENDENCIES ${DEPENDENCY})
    else()
      list(APPEND REQUIRED_DEPENDENCIES ${DEPENDENCY})
    endif()
  endforeach()

  # Declare dependencies dependent on some configuration settings
  if(IGRAPH_GLPK_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES GLPK)
  endif()
  if(IGRAPH_GMP_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES GMP)
  endif()
  if(IGRAPH_GRAPHML_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES LibXml2)
  endif()

  # Find dependencies
  foreach(DEPENDENCY ${REQUIRED_DEPENDENCIES} ${OPTIONAL_DEPENDENCIES})
    find_package(${DEPENDENCY})
  endforeach()

  # Export some aliases that will be used in config.h
  set(HAVE_GLPK ${GLPK_FOUND})
  set(HAVE_GMP ${GMP_FOUND})
  set(HAVE_LIBXML ${LIBXML2_FOUND})

  # Override libraries of vendored dependencies even if they were somehow
  # detected above
  foreach(DEPENDENCY ${VENDORABLE_DEPENDENCIES})
    if(IGRAPH_USE_INTERNAL_${DEPENDENCY})
      string(TOLOWER "${DEPENDENCY}" LIBNAME)
      set(${DEPENDENCY}_LIBRARIES ${LIBNAME}_vendored)
    endif()
  endforeach()

  message(STATUS ${BLAS_LIBRARIES})
endmacro()
