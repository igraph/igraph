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
)

find_library(PLFIT_LIBRARY
  NAMES plfit
)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  PLFIT
  DEFAULT_MSG
  PLFIT_LIBRARY
  PLFIT_INCLUDE_DIR
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(PLFIT_INCLUDE_DIR)
mark_as_advanced(PLFIT_LIBRARY)

if(PLFIT_FOUND)
  set(PLFIT_LIBRARIES ${PLFIT_LIBRARY})
  set(PLFIT_INCLUDE_DIRS ${PLFIT_INCLUDE_DIR})
endif()
