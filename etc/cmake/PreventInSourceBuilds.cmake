# Original source of this script:
# https://raw.githubusercontent.com/InsightSoftwareConsortium/ITK/master/CMake/PreventInSourceBuilds.cmake
#
# Thanks to the ITK project!
#
# This function will prevent in-source builds
function(AssureOutOfSourceBuilds)
  # make sure the user doesn't play dirty with symlinks
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if("${srcdir}" STREQUAL "${bindir}")
    message("##########################################################################")
    message("# igraph should not be configured & built in the igraph source directory")
    message("# You must run cmake in a build directory.")
	message("#")
	message("# Example:")
    message("# mkdir build; cd build; cmake ..; make")
    message("#")
    message("# NOTE: Given that you already tried to make an in-source build")
    message("#       CMake have already created several files & directories")
	message("#       in your source tree. If you are using git, run 'git clean -dfx'")
	message("#       to start from scratch. If you don't have git, remove")
	message("#       CMakeCache.txt and the CMakeFiles/ folder from the top of")
	message("#       the source tree.")
    message("#")
    message("##########################################################################")
    message("")
    message(FATAL_ERROR "Quitting configuration")
  endif()
endfunction()

AssureOutOfSourceBuilds()
