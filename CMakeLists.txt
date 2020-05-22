cmake_minimum_required(VERSION 3.0)

# Make use of ccache if it is present on the host system
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

# Declare the project, its version number and language
project(igraph VERSION 0.9.0 LANGUAGES C CXX)

# Set C and C++ standard version
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add etc/cmake to CMake's search path so we can put our private stuff there
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/etc/cmake)

# Add switches to use sanitizers if needed
include(sanitizers)

# Add version information
# TODO(ntamas): get this from Git or IGRAPH_VERSION
set(PACKAGE_VERSION "0.9.0-cmake-experimental")
set(PACKAGE_VERSION_MAJOR 0)
set(PACKAGE_VERSION_MINOR 9)
set(PACKAGE_VERSION_PATCH 0)
set(PACKAGE_VERSION_PRERELEASE "cmake-experimental")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/igraph_version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/igraph_version.h
)

# Run compile-time checks, generate config.h and igraph_threading.h
include(CheckFunctionExists)
include(CheckSymbolExists)
check_function_exists(stpcpy HAVE_STPCPY)
check_symbol_exists(stpcpy "string.h" HAVE_STPCPY_SIGNATURE)
check_symbol_exists(strdup "string.h" HAVE_STRDUP)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/src/config.h
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/igraph_threading.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/igraph_threading.h
)

# Find dependencies
find_package(ARPACK REQUIRED)
find_package(CXSparse REQUIRED)
find_package(LAPACK REQUIRED)

# Enable unit tests. Behave nicely and do this only if we are not being
# included as a sub-project in another CMake project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

# Traverse subdirectories
add_subdirectory(src)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(tests)
endif()

# Show result of configuration
message(STATUS " ")
message(STATUS "----------------------------------")
message(STATUS "igraph configured successfully.")
message(STATUS "----------------------------------")
message(STATUS " ")
