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

#ifndef REST_TYPES_H
#define REST_TYPES_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif

#include "igraph_error.h"
#include <stddef.h>
#include <math.h>
#include <stdio.h>

typedef int    igraph_integer_t;
typedef double igraph_real_t;
typedef int    igraph_bool_t;

/* Replacements for printf that print doubles in the same way on all platforms
 * (even for NaN and infinities) */
int igraph_real_printf(igraph_real_t val);
int igraph_real_fprintf(FILE *file, igraph_real_t val);
int igraph_real_snprintf(char* str, size_t size, igraph_real_t val);

/* Replacements for printf that print doubles in the same way on all platforms
 * (even for NaN and infinities) with the largest possible precision */
int igraph_real_printf_precise(igraph_real_t val);
int igraph_real_fprintf_precise(FILE *file, igraph_real_t val);
int igraph_real_snprintf_precise(char* str, size_t size, igraph_real_t val);

/* igraph_i_fdiv is needed here instead of in igraph_math.h because
 * some constants use it */
double igraph_i_fdiv(const double a, const double b);

#if defined(INFINITY)
#  define IGRAPH_INFINITY INFINITY
#  define IGRAPH_POSINFINITY INFINITY
#  define IGRAPH_NEGINFINITY (-INFINITY)
#else
#  define IGRAPH_INFINITY (igraph_i_fdiv(1.0, 0.0))
#  define IGRAPH_POSINFINITY (igraph_i_fdiv(1.0, 0.0))
#  define IGRAPH_NEGINFINITY (igraph_i_fdiv(-1.0, 0.0))
#endif

int igraph_finite(double x);
#define IGRAPH_FINITE(x) igraph_finite(x)

int igraph_is_nan(double x);
int igraph_is_inf(double x);
int igraph_is_posinf(double x);
int igraph_is_neginf(double x);

#if defined(NAN)
#  define IGRAPH_NAN NAN
#elif defined(INFINITY)
#  define IGRAPH_NAN (INFINITY/INFINITY)
#else
#  define IGRAPH_NAN (igraph_i_fdiv(0.0, 0.0))
#endif

__END_DECLS

#endif
