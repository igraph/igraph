# CMake script that generates the IGRAPH_VERSION file in the build folder
#
# Script variables that need to be set before calling it via "cmake -P":
#
# * IGRAPH_VERSION should be set to the exact version number
# * VERSION_FILE_PATH should be set to the name of the version file

FILE(WRITE "${VERSION_FILE_PATH}" "${IGRAPH_VERSION}")

