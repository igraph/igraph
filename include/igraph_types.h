/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_TYPES_H
#define IGRAPH_TYPES_H

#include "igraph_decls.h"

__BEGIN_DECLS

#ifdef __cplusplus
    #define __STDC_FORMAT_MACROS   /* needed for PRId32 and PRId64 from inttypes.h on Linux */
#endif

#include "igraph_config.h"

#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>


#if !defined(IGRAPH_INTEGER_SIZE)
#  error "igraph integer size not defined; check the value of IGRAPH_INTEGER_SIZE when compiling"
#elif IGRAPH_INTEGER_SIZE == 64
typedef int64_t igraph_integer_t;
typedef uint64_t igraph_uint_t;
#elif IGRAPH_INTEGER_SIZE == 32
typedef int32_t igraph_integer_t;
typedef uint32_t igraph_uint_t;
#else
#  error "Invalid igraph integer size; check the value of IGRAPH_INTEGER_SIZE when compiling"
#endif

typedef double igraph_real_t;

/* IGRAPH_BOOL_TYPE is set to 'bool' by default, and it is not meant to be
 * overridden, except for the R interface where we know what we are doing.
 * See igraph_config.h for more info */
typedef IGRAPH_BOOL_TYPE igraph_bool_t;

/* printf format specifier for igraph_integer_t */
#if IGRAPH_INTEGER_SIZE == 64
#  define IGRAPH_PRId PRId64
#  define IGRAPH_PRIu PRIu64
#else
#  define IGRAPH_PRId PRId32
#  define IGRAPH_PRIu PRIu32
#endif

/* maximum and minimum allowed values for igraph_integer_t */
#if IGRAPH_INTEGER_SIZE == 64
#  define IGRAPH_INTEGER_MAX INT64_MAX
#  define IGRAPH_INTEGER_MIN INT64_MIN
#else
#  define IGRAPH_INTEGER_MAX INT32_MAX
#  define IGRAPH_INTEGER_MIN INT32_MIN
#endif

/* maximum and minimum allowed values for igraph_uint_t */
#if IGRAPH_INTEGER_SIZE == 64
#  define IGRAPH_UINT_MAX UINT64_MAX
#  define IGRAPH_UINT_MIN UINT64_MIN
#else
#  define IGRAPH_UINT_MAX UINT32_MAX
#  define IGRAPH_UINT_MIN UINT32_MIN
#endif


/**
 * \define IGRAPH_VCOUNT_MAX
 * \brief The maximum number of vertices supported in igraph graphs.
 *
 * The value of this constant is one less than \c IGRAPH_INTEGER_MAX .
 * When igraph is compiled in 32-bit mode, this means that you are limited
 * to 2<superscript>31</superscript> – 2 (about 2.1 billion) vertices. In
 * 64-bit mode, the limit is 2<superscript>63</superscript> – 2 so you are much
 * more likely to hit out-of-memory issues due to other reasons before reaching
 * this limit.
 */
#define IGRAPH_VCOUNT_MAX (IGRAPH_INTEGER_MAX-1)
/* The 'os' and 'is' vectors in igraph_t have vcount+1 elements,
 * thus this cannot currently be larger than IGRAPH_INTEGER_MAX-1
 */

/**
 * \define IGRAPH_ECOUNT_MAX
 * \brief The maximum number of edges supported in igraph graphs.
 *
 * The value of this constant is half of \c IGRAPH_INTEGER_MAX .
 * When igraph is compiled in 32-bit mode, this means that you are limited
 * to approximately 2<superscript>30</superscript> (about 1.07 billion)
 * vertices. In 64-bit mode, the limit is approximately
 * 2<superscript>62</superscript> so you are much more likely to hit
 * out-of-memory issues due to other reasons before reaching this limit.
 */
#define IGRAPH_ECOUNT_MAX (IGRAPH_INTEGER_MAX/2)
/* The endpoints of edges are often stored in a vector twice the length
 * of the edge count, thus this cannot be larger than IGRAPH_INTEGER_MAX/2.
 * Some of the overflow checking code relies on this. */

/* Replacements for printf that print doubles in the same way on all platforms
 * (even for NaN and infinities) */
IGRAPH_EXPORT int igraph_real_printf(igraph_real_t val);
IGRAPH_EXPORT int igraph_real_fprintf(FILE *file, igraph_real_t val);
IGRAPH_EXPORT int igraph_real_printf_aligned(int width, igraph_real_t val);
IGRAPH_EXPORT int igraph_real_fprintf_aligned(FILE *file, int width, igraph_real_t val);
IGRAPH_EXPORT int igraph_real_snprintf(char *str, size_t size, igraph_real_t val);

/* Replacements for printf that print doubles in the same way on all platforms
 * (even for NaN and infinities) with the largest possible precision */
IGRAPH_EXPORT int igraph_real_printf_precise(igraph_real_t val);
IGRAPH_EXPORT int igraph_real_fprintf_precise(FILE *file, igraph_real_t val);
IGRAPH_EXPORT int igraph_real_snprintf_precise(char *str, size_t size, igraph_real_t val);

#define IGRAPH_INFINITY ((double)INFINITY)
#define IGRAPH_POSINFINITY IGRAPH_INFINITY
#define IGRAPH_NEGINFINITY (-IGRAPH_INFINITY)

IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_finite(double x);
#define IGRAPH_FINITE(x) igraph_finite(x)

IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_is_nan(double x);
IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_is_inf(double x);
IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_is_posinf(double x);
IGRAPH_DEPRECATED IGRAPH_EXPORT int igraph_is_neginf(double x);

#define IGRAPH_NAN ((double)NAN)

__END_DECLS

#endif
