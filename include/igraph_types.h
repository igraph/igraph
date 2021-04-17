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

#include "igraph_error.h"
#include <stddef.h>
#include <math.h>
#include <stdio.h>

typedef int    igraph_integer_t;
typedef double igraph_real_t;
typedef int    igraph_bool_t;

/* printf format specifier for igraph_integer_t */
#define IGRAPH_PRId "d"

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
