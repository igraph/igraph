/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2008  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#include "config.h"

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

/*
 * Compiler-related hacks, mostly because of Microsoft Visual C++
 */
double igraph_i_fdiv(const double a, const double b);
int igraph_i_snprintf(char *buffer, size_t count, const char *format, ...);
double igraph_i_round(double X);

double igraph_log2(const double a);
long double igraph_logbl(long double a);
double igraph_log1p(double a);
double igraph_fmin(double a, double b);
#ifndef HAVE_LOG2
#define log2(a) igraph_log2(a)
#endif
#ifndef HAVE_LOGBL
#define logbl(a) igraph_logbl(a)
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
#  define M_PI 3.14159265358979323846
#endif
#if !defined(M_LN2)
#  define M_LN2 0.69314718055994530942
#endif

__END_DECLS

#endif

