include(GetGitRevisionDescription)

# At this point, igraph is either the main CMake project or a subproject of
# another project. CMAKE_SOURCE_DIR would point to the root of the main
# project if we are a subproject so we cannot use that; we need to use
# CMAKE_CURRENT_SOURCE_DIR to get the directory containing the CMakeLists.txt
# file that version.cmake was included from, which is the top-level
# CMakeLists.txt file of igraph itself
set(VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/IGRAPH_VERSION")
set(NEXT_VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/NEXT_VERSION")

if(EXISTS "${VERSION_FILE}")
  file(READ "${VERSION_FILE}" PACKAGE_VERSION)
  string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
  message(STATUS "Version number: ${PACKAGE_VERSION}")
else()
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
  elseif(EXISTS "${NEXT_VERSION_FILE}")
    file(READ "${NEXT_VERSION_FILE}" PACKAGE_VERSION)
    string(STRIP "${PACKAGE_VERSION}" PACKAGE_VERSION)
    string(APPEND PACKAGE_VERSION "-dev")
    message(STATUS "Version number: ${PACKAGE_VERSION}")
  else()
    message(STATUS "Cannot find out the version number of this package; IGRAPH_VERSION is missing.")
    message(STATUS "")
    message(STATUS "The official igraph tarballs should contain this file, therefore you are")
    message(STATUS "most likely trying to compile a development version yourself. The development")
    message(STATUS "versions need Git to be able to determine the version number of igraph.")
    message(STATUS "")
    if(Git_FOUND)
      message(STATUS "It seems like you do have Git but it failed to determine the package version number.")
      message(STATUS "")
      message(STATUS "Git was found at: ${GIT_EXECUTABLE}")
      message(STATUS "The version number detection failed with: ${PACKAGE_VERSION}")
      message(STATUS "")
      message(STATUS "Most frequently this is caused by a shallow Git checkout that contains no tags in the history.")
    else()
      message(STATUS "Please install Git, make sure it is in your path, and then try again.")
    endif()
    message(STATUS "")
    message(FATAL_ERROR "Configuration failed.")
  endif()
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
# folder, for the sake of creating a tarball. This is needed only if igraph is
# the main project
if(NOT PROJECT_NAME)
  add_custom_target(
    versionfile
    BYPRODUCTS "${CMAKE_BINARY_DIR}/IGRAPH_VERSION"
    COMMAND "${CMAKE_COMMAND}"
      -DIGRAPH_VERSION="${PACKAGE_VERSION}"
      -DVERSION_FILE_PATH="${CMAKE_BINARY_DIR}/IGRAPH_VERSION"
      -P "${CMAKE_SOURCE_DIR}/etc/cmake/create_igraph_version_file.cmake"
    COMMENT "Generating IGRAPH_VERSION file in build folder"
  )
endif()
