# Declare the files needed to compile our vendored plfit copy
add_library(
  qhull_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  # list C files here
)

target_include_directories(
  qhull_vendored
  #PRIVATE
  #${PROJECT_SOURCE_DIR}/include
  #${PROJECT_BINARY_DIR}/include
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET qhull_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(qhull_vendored PRIVATE IGRAPH_STATIC)

use_all_warnings(qhull_vendored)
