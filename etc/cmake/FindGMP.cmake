# Inspired by http://code.google.com/p/origin/source/browse/trunk/cmake/FindGMP.cmake

# Copyright (c) 2008-2010 Kent State University
# Copyright (c) 2011-2012 Texas A&M University
#
# This file is distributed under the MIT License. See
# http://www.opensource.org/licenses/mit-license.php for terms and conditions.
#
# Some modifications made by Tamas Nepusz to ensure that the module fits better
# with the de facto conventions of FindXXX.cmake scripts

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
)

find_library(GMP_LIBRARY
  NAMES gmp
)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "GMP"
  DEFAULT_MSG
  GMP_LIBRARY
  GMP_INCLUDE_DIR
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(GMP_INCLUDE_DIR)
mark_as_advanced(GMP_LIBRARY)

if(GMP_FOUND)
  set(GMP_LIBRARIES ${GMP_LIBRARY})
  set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
endif()
