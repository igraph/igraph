include(PadString)

macro(find_dependencies)
  # Declare configuration options for dependencies
  option(IGRAPH_GLPK_SUPPORT "Compile igraph with GLPK support" YES)
  option(IGRAPH_GMP_SUPPORT "Compile igraph with GMP support" YES)
  option(IGRAPH_GRAPHML_SUPPORT "Compile igraph with GraphML support" YES)

  # Declare dependencies
  set(REQUIRED_DEPENDENCIES ARPACK CXSparse LAPACK)
  set(OPTIONAL_DEPENDENCIES FLEX BISON)

  # Declare dependencies dependent on some configuration settings
  if(IGRAPH_GLPK_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES GLPK)
  endif(IGRAPH_GLPK_SUPPORT)
  if(IGRAPH_GMP_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES GMP)
  endif(IGRAPH_GMP_SUPPORT)
  if(IGRAPH_GRAPHML_SUPPORT)
    list(APPEND REQUIRED_DEPENDENCIES LibXml2)
  endif(IGRAPH_GRAPHML_SUPPORT)

  # Find dependencies
  foreach(DEPENDENCY ${REQUIRED_DEPENDENCIES} ${OPTIONAL_DEPENDENCIES})
    find_package(${DEPENDENCY})
  endforeach()

  # Export some aliases that will be used in config.h
  set(HAVE_GLPK ${GLPK_FOUND})
  set(HAVE_GMP ${GMP_FOUND})
  set(HAVE_LIBXML ${LIBXML2_FOUND})
endmacro()
