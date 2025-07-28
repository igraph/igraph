file(GLOB_RECURSE infomap_src CONFIGURE_DEPENDS "src/*.cpp")

add_library(infomap_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  ${infomap_src}
)

target_include_directories(
  infomap_vendored

  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src/

  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_BINARY_DIR}/src # config.h
)

if(BUILD_SHARED_LIBS)
  set_property(TARGET infomap_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(infomap_vendored PRIVATE IGRAPH_STATIC)

# Infomap uses unsigned indices in OpenMP loops, which requires
# at least OpenMP 3.0. As of 2025, MSVC only supports OpenMP 2.0,
# thus this check is necessary.
if(IGRAPH_OPENMP_SUPPORT AND OpenMP_CXX_VERSION VERSION_GREATER_EQUAL 3)
  target_link_libraries(infomap_vendored PRIVATE OpenMP::OpenMP_CXX)
endif()
