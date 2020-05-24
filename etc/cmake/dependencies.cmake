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

function(print_bool VAR HEADING)
  if(${VAR})
    set(LABEL "yes")
  else()
    set(LABEL "no")
  endif()
  print_str(${VAR} ${HEADING} ${LABEL})
endfunction()

function(print_str VAR HEADING LABEL)
  string(LENGTH "${HEADING}" HEADING_LENGTH)
  math(EXPR REMAINING_WIDTH "30 - ${HEADING_LENGTH}")
  pad_string(PADDED ${REMAINING_WIDTH} " " "${LABEL}")
  message(STATUS "${HEADING}: ${PADDED}")
endfunction()

macro(summarize_dependencies)
  set(ALL_DEPENDENCIES ${REQUIRED_DEPENDENCIES} ${OPTIONAL_DEPENDENCIES})
  list(SORT ALL_DEPENDENCIES CASE INSENSITIVE)

  message(STATUS " ")
  message(STATUS "-----[ Build configuration ]----")
  print_str(PACKAGE_VERSION "Version" "${PACKAGE_VERSION}")
  print_str(CMAKE_BUILD_TYPE "CMake build type" "${CMAKE_BUILD_TYPE}")

  message(STATUS " ")
  message(STATUS "----------[ Features ]----------")
  print_bool(IGRAPH_GLPK_SUPPORT "GLPK for optimization")
  print_bool(IGRAPH_GMP_SUPPORT "GMP for big integers")
  print_bool(IGRAPH_GRAPHML_SUPPORT "Reading GraphML files")
  message(STATUS " ")

  message(STATUS "--------[ Dependencies ]--------")
  foreach(DEPENDENCY ${ALL_DEPENDENCIES})
    print_bool(${DEPENDENCY}_FOUND "${DEPENDENCY}")
  endforeach()
  message(STATUS " ")

  set(MISSING_DEPENDENCIES)
  foreach(DEPENDENCY ${REQUIRED_DEPENDENCIES})
    if(NOT ${DEPENDENCY}_FOUND)
      list(APPEND MISSING_DEPENDENCIES ${DEPENDENCY})
    endif()
  endforeach()

  if(MISSING_DEPENDENCIES)
    list(JOIN MISSING_DEPENDENCIES ", " GLUED)
    message(FATAL_ERROR "The following dependencies are missing: ${GLUED}")
  else()
    message(STATUS "igraph configured successfully.")
    message(STATUS " ")
  endif()
endmacro()
