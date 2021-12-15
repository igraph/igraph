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

#ifndef _GNU_SOURCE
    #define _GNU_SOURCE 1
#endif

#ifdef __cplusplus
    #define __STDC_FORMAT_MACROS   /* needed for PRId32 and PRId64 from inttypes.h on Linux */
#endif

#include "igraph_config.h"
#include "igraph_error.h"

#include <inttypes.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>

#if !defined(IGRAPH_INTEGER_SIZE)
#  error "igraph integer size not defined; check the value of IGRAPH_INTEGER_SIZE when compiling"
#elif IGRAPH_INTEGER_SIZE == 64
typedef int64_t igraph_integer_t;
#elif IGRAPH_INTEGER_SIZE == 32
typedef int32_t igraph_integer_t;
#else
#  error "Invalid igraph integer size; check the value of IGRAPH_INTEGER_SIZE when compiling"
#endif

typedef double igraph_real_t;
typedef int    igraph_bool_t;

/* printf format specifier for igraph_integer_t */
#if IGRAPH_INTEGER_SIZE == 64
#  define IGRAPH_PRId PRId64
#else
#  define IGRAPH_PRId PRId32
#endif

/* maximum allowed value for igraph_integer_t */
#if IGRAPH_INTEGER_SIZE == 64
#  define IGRAPH_INTEGER_MAX INT64_MAX
#else
#  define IGRAPH_INTEGER_MAX INT32_MAX
#endif

/* Replacements for printf that print doubles in the same way on all platforms
 * (even for NaN and infinities) */
IGRAPH_EXPORT int igraph_real_printf(igraph_real_t val);
IGRAPH_EXPORT int igraph_real_fprintf(FILE *file, igraph_real_t val);
IGRAPH_EXPORT int igraph_real_snprintf(char* str, size_t size, igraph_real_t val);

/* Replacements for printf that print doubles in the same way on all platforms
 * (even for NaN and infinities) with the largest possible precision */
IGRAPH_EXPORT int igraph_real_printf_precise(igraph_real_t val);
IGRAPH_EXPORT int igraph_real_fprintf_precise(FILE *file, igraph_real_t val);
IGRAPH_EXPORT int igraph_real_snprintf_precise(char* str, size_t size, igraph_real_t val);

#define IGRAPH_INFINITY INFINITY
#define IGRAPH_POSINFINITY INFINITY
#define IGRAPH_NEGINFINITY (-INFINITY)

IGRAPH_EXPORT int igraph_finite(double x);
#define IGRAPH_FINITE(x) igraph_finite(x)

IGRAPH_EXPORT int igraph_is_nan(double x);
IGRAPH_EXPORT int igraph_is_inf(double x);
IGRAPH_EXPORT int igraph_is_posinf(double x);
IGRAPH_EXPORT int igraph_is_neginf(double x);

#define IGRAPH_NAN NAN

__END_DECLS

#endif
