# https://raw.githubusercontent.com/dune-project/dune-istl/master/cmake/modules/FindARPACK.cmake
#
# This file is taken from:
#
# DUNE, the Distributed and Unified Numerics Environment
# GPLv2 licensed
#
# .. cmake_module::
#
#    Module that checks whether ARPACK is available and usable.
#
#    Variables used by this module which you may want to set:
#
#    :ref:`ARPACK_ROOT`
#       Path list to search for ARPACK.
#
#    Sets the following variables:
#
#    :code:`ARPACK_FOUND`
#       True if ARPACK available.
#
#    :code:`ARPACK_LIBRARIES`
#       Link against these libraries to use ARPACK.
#
# .. cmake_variable:: ARPACK_ROOT
#
#    You may set this variable to have :ref:`FindARPACK` look
#    for the ARPACK package in the given path before inspecting
#    system paths.
#

# look for library, only at positions given by the user
find_library(ARPACK_LIBRARY
  NAMES "arpack"
  PATHS ${ARPACK_PREFIX} ${ARPACK_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  NO_DEFAULT_PATH
)

# look for library files, including default paths
find_library(ARPACK_LIBRARY
  NAMES "arpack"
  PATH_SUFFIXES "lib" "lib32" "lib64"
)

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND if the
# searches above were not successful; without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they
# are set to NOTFOUND. Please set them or make sure they are set and tested
# correctly in the CMake files."
if(ARPACK_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ARPACK_LIBRARY})
endif()

# end of header usability check
cmake_pop_check_state()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ARPACK"
  DEFAULT_MSG
  ARPACK_LIBRARY
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(ARPACK_LIBRARY)

# if headers are found, store results
if(ARPACK_FOUND)
  set(ARPACK_LIBRARIES ${ARPACK_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ARPACK succeeded:\n"
    "Libraries to link against: ${ARPACK_LIBRARIES}\n\n")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of ARPACK failed:\n"
    "Libraries to link against: ${ARPACK_LIBRARIES}\n\n")
endif()
