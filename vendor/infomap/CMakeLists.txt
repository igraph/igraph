# Declare the files needed to compile our vendored Infomap copy
add_library(infomap_vendored
  OBJECT
  EXCLUDE_FROM_ALL
  src/core/BiasedMapEquation.cpp
  src/core/InfoEdge.cpp
  src/core/InfoNode.cpp
  src/core/InfomapBase.cpp
  src/core/LossyMapEquation.cpp
  src/core/MemMapEquation.cpp
  src/core/MetaMapEquation.cpp
  src/core/RegularizedMultilayerMapEquation.cpp
  src/core/StateNetwork.cpp
  src/core/iterators/InfomapIterator.cpp
  src/io/ClusterMap.cpp
  src/io/Config.cpp
  src/io/Network.cpp
  src/io/NetworkBuilder.cpp
  src/io/Output.cpp
  src/io/OutputFormats.cpp
  src/io/OutputPlan.cpp
  src/io/OutputView.cpp
  src/io/ParameterCatalog.cpp
  src/io/ProgramInterface.cpp
  src/io/RunMetadata.cpp
  src/io/RunReport.cpp
  src/io/TrialResults.cpp
  src/utils/Console.cpp
  src/utils/FileURI.cpp
  src/utils/FlowCalculator.cpp
  src/utils/Log.cpp
  src/utils/MemoryUsage.cpp
  src/utils/ThreadConfig.cpp
  src/utils/convert.cpp
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

# Header-only dependencies of Infomap. fmt is reached from Infomap's own
# headers, so it must be visible wherever Infomap.h is included.
target_include_directories(
  infomap_vendored SYSTEM

  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/fmt/include

  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/nlohmann_json/include
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
