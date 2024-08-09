# Minimum CMake that we require is 3.18.
# Some of the recent features we use:
#  * --ignore-eol when comparing unit test results with expected outcomes (3.14)
#  * CROSSCOMPILING_EMULATOR can be a semicolon-separated list to pass arguments (3.15)
#  * SKIP_REGULAR_EXPRESSION to handle skipped tests properly (3.16)
#  * CheckLinkerFlag for HAVE_NEW_DTAGS test (3.18)
#  * cmake -E cat (3.18)
cmake_minimum_required(VERSION 3.18...3.27)

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
include(UseCCacheWhenInstalled)

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

# Detect is certain attributes are supported by the compiler
include(attribute_support)

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

# Enable fuzzer instrumentation if needed
# FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION is a conventional
# macro used to adapt code for fuzzability, for example by
# reducing largest allowed graph sizes when reading various
# file formats.
if(BUILD_FUZZING)
  add_compile_options(-fsanitize=fuzzer-no-link)
  add_compile_definitions(FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION)
endif()

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
include(CheckIncludeFiles)
include(CMakePushCheckState)

# First we check for some functions and symbols
cmake_push_check_state()
if(NEED_LINKING_AGAINST_LIBM)
  list(APPEND CMAKE_REQUIRED_LIBRARIES m)
endif()
check_symbol_exists(strcasecmp strings.h HAVE_STRCASECMP)
check_symbol_exists(strncasecmp strings.h HAVE_STRNCASECMP)
check_symbol_exists(_stricmp string.h HAVE__STRICMP)
check_symbol_exists(_strnicmp string.h HAVE__STRNICMP)
check_symbol_exists(strdup string.h HAVE_STRDUP)
check_symbol_exists(strndup string.h HAVE_STRNDUP)
check_include_files(xlocale.h HAVE_XLOCALE)
if(HAVE_XLOCALE)
  # On BSD, uselocale() is in xlocale.h instead of locale.h.
  # Some systems provide xlocale.h, but uselocale() is still in locale.h,
  # thus we try both.
  check_symbol_exists(uselocale "xlocale.h;locale.h" HAVE_USELOCALE)
else()
  check_symbol_exists(uselocale locale.h HAVE_USELOCALE)
endif()
check_symbol_exists(_configthreadlocale locale.h HAVE__CONFIGTHREADLOCALE)
cmake_pop_check_state()

# Check for 128-bit integer multiplication support, floating-point endianness,
# support for built-in overflow detection and fast bit operation support.
include(ieee754_endianness)
include(uint128_support)
include(bit_operations_support)
include(safe_math_support)

if(NOT HAVE_USELOCALE AND NOT HAVE__CONFIGTHREADLOCALE)
  message(WARNING "igraph cannot set per-thread locale on this platform. igraph_enter_safelocale() and igraph_exit_safelocale() will not be safe to use in multithreaded programs.")
endif()

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
    # directory that somehow gets into the coverage results. /Applications and
    # /Library/Developer are for macOS -- they exclude files from the macOS SDK.
    EXCLUDE "io/*.l" "src/io/parsers/*" "io/parsers/*" "/Applications/Xcode*" "/Library/Developer/*" "examples/*" "interfaces/*" "tests/*" "vendor/pcg/*"
  )
endif()

# Generate configuration headers
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/src/config.h
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/igraph_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/igraph_config.h
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

# Traverse subdirectories. vendor/ should come first because code in
# src/CMakeLists.txt depends on targets in vendor/
add_subdirectory(vendor)
add_subdirectory(src)
add_subdirectory(interfaces)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(tests)
endif()
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_FUZZING)
  add_subdirectory(fuzzing)
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
