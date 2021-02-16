include(GetGitRevisionDescription)

set(VERSION_FILE "${CMAKE_SOURCE_DIR}/IGRAPH_VERSION")
set(NEXT_VERSION_FILE "${CMAKE_SOURCE_DIR}/NEXT_VERSION")

find_package(Git QUIET)
if(Git_FOUND)
  git_describe(PACKAGE_VERSION)
else()
  set(PACKAGE_VERSION "NOTFOUND")
endif()

if(PACKAGE_VERSION)
  if(EXISTS "${NEXT_VERSION_FILE}")
    file(READ "${NEXT_VERSION_FILE}" PACKAGE_VERSION)
    string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
    get_git_head_revision(GIT_REFSPEC GIT_COMMIT_HASH)
    string(SUBSTRING "${GIT_COMMIT_HASH}" 0 8 GIT_COMMIT_HASH_SHORT)
  string(APPEND PACKAGE_VERSION "-dev+${GIT_COMMIT_HASH_SHORT}")
  endif()
  message(STATUS "Version number from Git: ${PACKAGE_VERSION}")
elseif(EXISTS "${VERSION_FILE}")
  file(READ "${VERSION_FILE}" PACKAGE_VERSION)
  string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
  message(STATUS "Version number: ${PACKAGE_VERSION}")
elseif(EXISTS "${NEXT_VERSION_FILE}")
  file(READ "${NEXT_VERSION_FILE}" PACKAGE_VERSION)
  string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
  string(APPEND PACKAGE_VERSION "-dev")
  message(STATUS "Version number: ${PACKAGE_VERSION}")
elseif(DEFINED ENV{TRAVIS_COMMIT})
  # git_describe() might fail on Travis if the latest git tag is too far away
  # in the past and Travis did not fetch it as part of the commit history. We
  # can fall back to using TRAVIS_COMMIT, though
  set(GIT_COMMIT_HASH "${ENV{TRAVIS_COMMIT}}")
  string(SUBSTRING "${GIT_COMMIT_HASH}" 0 8 GIT_COMMIT_HASH_SHORT)
  if(EXISTS "${NEXT_VERSION_FILE}")
    file(READ "${NEXT_VERSION_FILE}" PACKAGE_VERSION)
    string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
  string(APPEND PACKAGE_VERSION "-dev+${GIT_COMMIT_HASH_SHORT}")
  elseif(EXISTS "${VERSION_FILE}")
    file(READ "${VERSION_FILE}" PACKAGE_VERSION)
    string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
  string(APPEND PACKAGE_VERSION "-post+${GIT_COMMIT_HASH_SHORT}")
  else()
    message(STATUS "Cannot find out the version number of this package; IGRAPH_VERSION and NEXT_VERSION are both missing.")
    message(STATUS "")
    message(FATAL_ERROR "Configuration failed.")
  endif()
  message(STATUS "Version number from Travis environment: ${PACKAGE_VERSION}")
else()
  message(STATUS "Cannot find out the version number of this package; IGRAPH_VERSION is missing.")
  message(STATUS "")
  message(STATUS "The official igraph tarballs should contain this file, therefore you are")
  message(STATUS "most likely trying to compile a development version yourself. The development")
  message(STATUS "versions need Git to be able to determine the version number of igraph.")
  message(STATUS "")
  message(STATUS "Please install Git, make sure it is in your path, and then try again.")
  message(STATUS "")
  message(FATAL_ERROR "Configuration failed.")
endif()

string(REGEX MATCH "^[^-]+" PACKAGE_VERSION_BASE "${PACKAGE_VERSION}")
string(
  REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9+])" "\\1;\\2;\\3"
  PACKAGE_VERSION_PARTS "${PACKAGE_VERSION_BASE}"
)
list(GET PACKAGE_VERSION_PARTS 0 PACKAGE_VERSION_MAJOR)
list(GET PACKAGE_VERSION_PARTS 1 PACKAGE_VERSION_MINOR)
list(GET PACKAGE_VERSION_PARTS 2 PACKAGE_VERSION_PATCH)

if(PACKAGE_VERSION MATCHES "^[^-]+-")
  string(
    REGEX REPLACE "^[^-]+-([^+]*)" "\\1" PACKAGE_VERSION_PRERELEASE "${PACKAGE_VERSION}"
  )
else()
  set(PACKAGE_VERSION_PRERELEASE "cmake-experimental")
endif()

# Add a target that we can use to generate an IGRAPH_VERSION file in the build
# folder, for the sake of creating a tarball
add_custom_target(
  versionfile
  BYPRODUCTS "${CMAKE_BINARY_DIR}/IGRAPH_VERSION"
  COMMAND "${CMAKE_COMMAND}"
    -DIGRAPH_VERSION="${PACKAGE_VERSION}"
  -DVERSION_FILE_PATH="${CMAKE_BINARY_DIR}/IGRAPH_VERSION"
  -P "${CMAKE_SOURCE_DIR}/etc/cmake/create_igraph_version_file.cmake"
  COMMENT "Generating IGRAPH_VERSION file in build folder"
)
