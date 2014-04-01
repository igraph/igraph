/* ========================================================================== */
/* === Include/cholmod_io64 ================================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_io64.h.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod_io64.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Definitions required for large file I/O, which must come before any other
 * #includes.  These are not used if -DNLARGEFILE is defined at compile time.
 * Large file support may not be portable across all platforms and compilers;
 * if you encounter an error here, compile your code with -DNLARGEFILE.  In
 * particular, you must use -DNLARGEFILE for MATLAB 6.5 or earlier (which does
 * not have the io64.h include file).
 */

#ifndef CHOLMOD_IO_H
#define CHOLMOD_IO_H

/* skip all of this if NLARGEFILE is defined at the compiler command line */
#ifndef NLARGEFILE

#if defined(MATLAB_MEX_FILE) || defined(MATHWORKS)

/* CHOLMOD is being compiled as a MATLAB mexFunction, or for use in MATLAB */
#include "io64.h"

#else

/* CHOLMOD is being compiled in a stand-alone library */
#undef  _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#undef  _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64

#endif

#endif

#endif

