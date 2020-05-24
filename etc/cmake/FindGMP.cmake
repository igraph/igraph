# Inspired by http://code.google.com/p/origin/source/browse/trunk/cmake/FindGMP.cmake

# Copyright (c) 2008-2010 Kent State University
# Copyright (c) 2011-2012 Texas A&M University
#
# This file is distributed under the MIT License. See
# http://www.opensource.org/licenses/mit-license.php for terms and conditions.

set(GMP_PREFIX "" CACHE PATH "Path to GMP prefix")

find_path(GMP_INCLUDE_DIR gmp.h PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include)
find_library(GMP_LIBRARY libgmp.a PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
  set(GMP_FOUND TRUE)
endif()

if(GMP_FOUND)
  MESSAGE(STATUS "Found GMP: ${GMP_INCLUDE_DIR}/gmp.h and ${GMP_LIBRARY}")
else()
  MESSAGE(STATUS "GMP not found")
endif()
