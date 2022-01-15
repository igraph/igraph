# Inspired by http://code.google.com/p/origin/source/browse/trunk/cmake/FindGMP.cmake

# Copyright (c) 2021 Tamas Nepusz
#
# This file is distributed under the MIT License. See
# http://www.opensource.org/licenses/mit-license.php for terms and conditions.
#
# Some modifications made by Tamas Nepusz to ensure that the module fits better
# with the de facto conventions of FindXXX.cmake scripts

find_path(PLFIT_INCLUDE_DIR
  NAMES plfit.h
  PATH_SUFFIXES plfit
)

find_library(PLFIT_LIBRARY
  NAMES plfit
)

# parse version from header
if(PLFIT_INCLUDE_DIR)
  set(PLFIT_VERSION_FILE ${PLFIT_INCLUDE_DIR}/plfit_version.h)
  file(READ ${PLFIT_VERSION_FILE} PLFIT_VERSION_FILE_CONTENTS)

  string(REGEX MATCH "#define[ ]+PLFIT_VERSION_MAJOR[ ]+[0-9]+"
    PLFIT_VERSION_MAJOR "${PLFIT_VERSION_FILE_CONTENTS}")
  string(REGEX REPLACE "#define[ ]+PLFIT_VERSION_MAJOR[ ]+([0-9]+)" "\\1"
    PLFIT_VERSION_MAJOR "${PLFIT_VERSION_MAJOR}")

  string(REGEX MATCH "#define[ ]+PLFIT_VERSION_MINOR[ ]+[0-9]+"
    PLFIT_VERSION_MINOR "${PLFIT_VERSION_FILE_CONTENTS}")
  string(REGEX REPLACE "#define[ ]+PLFIT_VERSION_MINOR[ ]+([0-9]+)" "\\1"
    PLFIT_VERSION_MINOR "${PLFIT_VERSION_MINOR}")

  string(REGEX MATCH "#define[ ]+PLFIT_VERSION_PATCH[ ]+[0-9]+"
    PLFIT_VERSION_PATCH "${PLFIT_VERSION_FILE_CONTENTS}")
  string(REGEX REPLACE "#define[ ]+PLFIT_VERSION_PATCH[ ]+([0-9]+)" "\\1"
    PLFIT_VERSION_PATCH "${PLFIT_VERSION_PATCH}")

  set(PLFIT_VERSION "${PLFIT_VERSION_MAJOR}.${PLFIT_VERSION_MINOR}.${PLFIT_VERSION_PATCH}")

  # compatibility variables
  set(PLFIT_VERSION_STRING "${PLFIT_VERSION}")
endif()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PLFIT
  FOUND_VAR PLFIT_FOUND
  REQUIRED_VARS
    PLFIT_LIBRARY
    PLFIT_INCLUDE_DIR
  VERSION_VAR PLFIT_VERSION
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(PLFIT_INCLUDE_DIR)
mark_as_advanced(PLFIT_LIBRARY)

if(PLFIT_FOUND)
  set(PLFIT_LIBRARIES ${PLFIT_LIBRARY})
  set(PLFIT_INCLUDE_DIRS ${PLFIT_INCLUDE_DIR})
endif()
