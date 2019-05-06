#ifndef BLISS_DEFS_HH
#define BLISS_DEFS_HH

#include <cassert>
#include <cstdarg>

#include "config.h"

/*
  Copyright (c) 2003-2015 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.
  
  This file is part of bliss.
  
  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_GMP == 1
#  define BLISS_USE_GMP
#endif

#ifdef USING_R
#include <R.h>
#define fatal_error(...) (error(__VA_ARGS__))
#endif

namespace bliss {

/**
 * The version number of bliss.
 */
static const char * const version = "0.73";

/*
 * If a fatal error (out of memory, internal error) is encountered,
 * this function is called.
 * There should not be a return from this function but exit or
 * a jump to code that deallocates the AbstractGraph instance that called this.
 */
#ifndef USING_R
void fatal_error(const char* fmt, ...);
#endif


#if defined(BLISS_DEBUG)
#define BLISS_CONSISTENCY_CHECKS
#define BLISS_EXPENSIVE_CONSISTENCY_CHECKS
#endif


#if defined(BLISS_CONSISTENCY_CHECKS)
/* Force a check that the found automorphisms are valid */
#define BLISS_VERIFY_AUTOMORPHISMS
#endif


#if defined(BLISS_CONSISTENCY_CHECKS)
/* Force a check that the generated partitions are equitable */
#define BLISS_VERIFY_EQUITABLEDNESS
#endif

} // namespace bliss



/*! \mainpage Bliss
 *
 * \section intro_sec Introduction
 *
 * This is the source code documentation of bliss,
 * produced by running <A href="http://www.doxygen.org">doxygen</A> in
 * the source directory.
 * The algorithms and data structures used in bliss are documented in
 * the papers found at the
 * <A href="http://www.tcs.hut.fi/Software/bliss">bliss web site</A>.
 *
 *
 * \section compile_sec Compiling
 *
 * Compiling bliss in Linux should be easy, just execute
 * \code
 * make
 * \endcode
 * in the bliss source directory.
 * This will produce the executable program \c bliss as well as
 * the library file \c libbliss.a that can be linked in other programs.
 * If you have the <A href="http://gmplib.org/">GNU Multiple Precision
 * Arithmetic Library</A> (GMP) installed in your machine, you can also use
 * \code
 * make gmp
 * \endcode
 * to enable exact computation of automorphism group sizes.
 *
 * When linking the bliss library \c libbliss.a in other programs,
 * remember to include the standard c++ library
 * (and the GMP library if you compiled bliss to include it).
 * For instance,
 * \code gcc -o test test.c -lstdc++ -lgmp -lbliss\endcode
 *
 * \section cppapi_sec The C++ language API
 *
 * The C++ language API is the main API to bliss;
 * all other APIs are just more or less complete variants of it.
 * The C++ API consists basically of the public methods in
 * the classes bliss::AbstractGraph, bliss::Graph, and bliss::Digraph.
 * For an example of its use,
 * see the \ref executable "source of the bliss executable".
 *
 *
 * \section capi_sec The C language API
 *
 * The C language API is given in the file bliss_C.h.
 * It is currently more restricted than the C++ API so
 * consider using the C++ API whenever possible.
 */


#endif
