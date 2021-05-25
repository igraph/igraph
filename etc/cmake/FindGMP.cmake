# Inspired by http://code.google.com/p/origin/source/browse/trunk/cmake/FindGMP.cmake

# Copyright (c) 2008-2010 Kent State University
# Copyright (c) 2011-2012 Texas A&M University
#
# This file is distributed under the MIT License. See
# http://www.opensource.org/licenses/mit-license.php for terms and conditions.
#
# Some modifications made by Tamas Nepusz to ensure that the module fits better
# with the de facto conventions of FindXXX.cmake scripts

set(GMP_PREFIX "" CACHE PATH "Path to GMP prefix")

find_path(GMP_INCLUDE_DIR gmp.h PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include)
find_library(GMP_LIBRARY gmp PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "GMP"
  DEFAULT_MSG
  GMP_LIBRARY
  GMP_INCLUDE_DIR
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(GMP_PREFIX)
mark_as_advanced(GMP_INCLUDE_DIR)
mark_as_advanced(GMP_LIBRARY)

if(GMP_FOUND)
  set(GMP_LIBRARIES ${GMP_LIBRARY})
endif()
