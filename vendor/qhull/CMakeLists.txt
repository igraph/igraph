# Declare the files needed to compile our vendored plfit copy
add_library(
  qhull_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  libqhull_r/accessors_r.c
  libqhull_r/geom_r.c
  libqhull_r/geom2_r.c
  libqhull_r/global_r.c
  libqhull_r/io_r.c
  libqhull_r/libqhull_r.c
  libqhull_r/mem_r.c
  libqhull_r/merge_r.c
  libqhull_r/poly_r.c
  libqhull_r/poly2_r.c
  libqhull_r/qset_r.c
  libqhull_r/random_r.c
  libqhull_r/rboxlib_r.c
  libqhull_r/stat_r.c
  libqhull_r/user_r.c
  libqhull_r/usermem_r.c
  libqhull_r/userprintf_r.c
  libqhull_r/userprintf_rbox_r.c

  libqhullcpp/Coordinates.cpp
  libqhullcpp/PointCoordinates.cpp
  libqhullcpp/Qhull.cpp
  libqhullcpp/QhullFacet.cpp
  libqhullcpp/QhullFacetList.cpp
  libqhullcpp/QhullFacetSet.cpp
  libqhullcpp/QhullHyperplane.cpp
  libqhullcpp/QhullPoint.cpp
  libqhullcpp/QhullPoints.cpp
  libqhullcpp/QhullPointSet.cpp
  libqhullcpp/QhullQh.cpp
  libqhullcpp/QhullRidge.cpp
  libqhullcpp/QhullSet.cpp
  libqhullcpp/QhullStat.cpp
  libqhullcpp/QhullUser.cpp
  libqhullcpp/QhullVertex.cpp
  libqhullcpp/QhullVertexSet.cpp
  libqhullcpp/qt-qhull.cpp
  libqhullcpp/RboxPoints.cpp
  libqhullcpp/RoadError.cpp
  libqhullcpp/RoadLogEvent.cpp
  libqhullcpp/usermem_r-cpp.cpp

)

target_include_directories(
  qhull_vendored
  PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include
  #PUBLIC
  #${CMAKE_CURRENT_SOURCE_DIR}
)

if (BUILD_SHARED_LIBS)
  set_property(TARGET qhull_vendored PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

# Since these are included as object files, they should call the
# function as is (without visibility specification)
target_compile_definitions(qhull_vendored PRIVATE IGRAPH_STATIC)

use_all_warnings(qhull_vendored)
