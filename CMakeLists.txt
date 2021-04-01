# Minimum CMake that we require is 3.16 because we use --ignore-eol when
# comparing unit test results with expected outcomes (added in 3.14) and we
# also use SKIP_REGULAR_EXPRESSION to handle skipped tests properly
cmake_minimum_required(VERSION 3.16)

# Add etc/cmake to CMake's search path so we can put our private stuff there
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/etc/cmake)

# Set a default build type if none was specified
# This must precede the project() line, which would set the CMAKE_BUILD_TYPE
# to 'Debug' with single-config generators on Windows.
# Note that we must do this only if PROJECT_NAME is not set at this point. If
# it is set, it means that igraph is being used as a subproject of another
# project.
if(NOT PROJECT_NAME)
  include(BuildType)
endif()

# Prevent in-source builds
include(PreventInSourceBuilds)

# Make use of ccache if it is present on the host system -- unless explicitly
# asked to disable it
option(USE_CCACHE "Use ccache to speed up compilation if it is installed" ON)
if(USE_CCACHE)
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  endif()
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

# Include some compiler-related helpers and set global compiler options
include(compilers)

# Set default symbol visibility to hidden
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Set C and C++ standard version
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
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

# Create configuration options for optional features
include(features)

# Handle dependencies and dependency-related configuration options
include(dependencies)
find_dependencies()

# Run compile-time checks, generate config.h and igraph_threading.h
include(CheckSymbolExists)

# The math library is needed on Linux to find symbols like exmp1 or fmin
if(UNIX)
  set(CMAKE_REQUIRED_LIBRARIES "m")
endif()

# First we check for some functions and symbols
check_symbol_exists(expm1 math.h HAVE_EXPM1)
check_symbol_exists(fmin math.h HAVE_FMIN)
check_symbol_exists(finite math.h HAVE_FINITE)
check_symbol_exists(isfinite math.h HAVE_ISFINITE)
check_symbol_exists(log2 math.h HAVE_LOG2)
check_symbol_exists(log1p math.h HAVE_LOG1P)
check_symbol_exists(rint math.h HAVE_RINT)
check_symbol_exists(rintf math.h HAVE_RINTF)
check_symbol_exists(round math.h HAVE_ROUND)
check_symbol_exists(stpcpy string.h HAVE_STPCPY)
check_symbol_exists(strcasecmp strings.h HAVE_STRCASECMP)
check_symbol_exists(strdup string.h HAVE_STRDUP)
check_symbol_exists(_stricmp string.h HAVE__STRICMP)

# Check for code coverage support
option(IGRAPH_ENABLE_CODE_COVERAGE "Enable code coverage calculation" OFF)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND IGRAPH_ENABLE_CODE_COVERAGE)
  include(CodeCoverage)
  append_coverage_compiler_flags()
  setup_target_for_coverage_lcov(
    NAME coverage
    EXECUTABLE "${CMAKE_COMMAND}" "--build" "${PROJECT_BINARY_DIR}" "--target" "check"
    # Generated files are excluded; apparently the CodeCoverage script has some
    # problems with them. Yes, the exclusion is correct, it refers to a nonexistent
    # directory that somehow gets into the coverage resolts. /Applications is for
    # macOS -- it excludes files from the macOS SDK.
    EXCLUDE "io/*.l" "io/parsers/*" "/Applications/Xcode*" "examples/*" "tests/*"
  )
endif()

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
  configure_file(
    ${PROJECT_SOURCE_DIR}/etc/cmake/CTestCustom.cmake.in
    ${PROJECT_BINARY_DIR}/CTestCustom.cmake
  )
endif()

# Traverse subdirectories
add_subdirectory(src)
add_subdirectory(vendor)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(doc)
endif()

# Configure packaging -- only if igraph is the top-level project and not a
# subproject
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(packaging)
endif()

# Show result of configuration
include(summary)
