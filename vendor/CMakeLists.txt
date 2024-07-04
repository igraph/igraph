add_subdirectory(cs)
add_subdirectory(f2c)
add_subdirectory(glpk)
add_subdirectory(lapack)
add_subdirectory(mini-gmp)
add_subdirectory(pcg)
add_subdirectory(plfit)

# We add the infomap vendored target here, so that
# we can keep infomap as a clean submodule.
add_library(infomap_vendored

  OBJECT
  EXCLUDE_FROM_ALL

  infomap/src/core/BiasedMapEquation.cpp
  infomap/src/core/InfoEdge.cpp
  infomap/src/core/InfomapBase.cpp
  infomap/src/core/InfoNode.cpp
  infomap/src/core/MemMapEquation.cpp
  infomap/src/core/MetaMapEquation.cpp
  infomap/src/core/StateNetwork.cpp

)

target_include_directories(
  infomap_vendored

  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}

  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/src # config.h for TLS
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET infomap_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(infomap_vendored PRIVATE IGRAPH_STATIC)
