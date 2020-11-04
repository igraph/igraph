# Minimum CMake that we require is 3.16 because we use --ignore-eol when
# comparing unit test results with expected outcomes (added in 3.14) and we
# also use SKIP_REGULAR_EXPRESSION to handle skipped tests properly
cmake_minimum_required(VERSION 3.16)

# Add etc/cmake to CMake's search path so we can put our private stuff there
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/etc/cmake)

# Make use of ccache if it is present on the host system
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

# Figure out the version number from Git
include(version)

# Declare the project, its version number and language
project(
  igraph
  VERSION ${PACKAGE_VERSION_BASE}
  DESCRIPTION "A library for creating and manipulating graphs"
  HOMEPAGE_URL https://igraph.org
  LANGUAGES C CXX
)

# Set C and C++ standard version
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Expose the BUILD_SHARED_LIBS option in the ccmake UI
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Add switches to use sanitizers and debugging helpers if needed
include(debugging)
include(sanitizers)

# Add version information
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/igraph_version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/igraph_version.h
)

# Handle dependencies and dependency-related configuration options
include(dependencies)
find_dependencies()

# Run compile-time checks, generate config.h and igraph_threading.h
include(CheckFunctionExists)
include(CheckSymbolExists)

# First we check for some functions and symbols
check_function_exists(stpcpy HAVE_STPCPY)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(_stricmp HAVE__STRICMP)
check_function_exists(expm1 HAVE_EXPM1)
check_function_exists(rint HAVE_RINT)
check_function_exists(rintf HAVE_RINTF)
check_symbol_exists(stpcpy "string.h" HAVE_STPCPY_SIGNATURE)
check_symbol_exists(strdup "string.h" HAVE_STRDUP)

# Check for TLS support and enable thread-local storage if possible
include(tls)

# Check for link-time optimization
include(lto)

# Generate configuration headers
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/src/config.h
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/igraph_threading.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/igraph_threading.h
)

# Enable unit tests. Behave nicely and do this only if we are not being
# included as a sub-project in another CMake project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

# Traverse subdirectories
add_subdirectory(src)
add_subdirectory(optional)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
add_subdirectory(doc)

# Configure packaging
include(packaging)

# Show result of configuration
include(summary)
