/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_MATH_H
#define IGRAPH_MATH_H

#include "igraph_decls.h"

#include "config.h"

#include <math.h>
#include <stddef.h>

__BEGIN_DECLS

/**
 * \def IGRAPH_SHORTEST_PATH_EPSILON
 *
 * Relative error threshold used in weighted shortest path calculations
 * to decide whether two shortest paths are of equal length.
 */
#define IGRAPH_SHORTEST_PATH_EPSILON 1e-10

/*
 * Compiler-related hacks, mostly because of Microsoft Visual C++
 */
double igraph_i_round(double X);
int igraph_i_snprintf(char *buffer, size_t count, const char *format, ...);

double igraph_log2(const double a);
double igraph_log1p(double a);
double igraph_fmin(double a, double b);
#ifndef HAVE_LOG2
    #define log2(a) igraph_log2(a)
#endif
#ifndef HAVE_LOG1P
    #define log1p(a) igraph_log1p(a)
#endif
#ifndef HAVE_FMIN
    #define fmin(a,b) igraph_fmin((a),(b))
#endif
#ifndef HAVE_ROUND
    #define round igraph_i_round
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif
#ifndef M_LN2
    #define M_LN2 0.69314718055994530942
#endif
#ifndef M_SQRT2
    #define M_SQRT2 1.4142135623730950488016887
#endif
#ifndef M_LN_SQRT_2PI
    #define M_LN_SQRT_2PI   0.918938533204672741780329736406 /* log(sqrt(2*pi))
    == log(2*pi)/2 */
#endif

int igraph_almost_equals(double a, double b, double eps);
int igraph_cmp_epsilon(double a, double b, double eps);

__END_DECLS

#endif
