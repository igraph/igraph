/* platform.h
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __PLATFORM_H__
#define __PLATFORM_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include <float.h>

__BEGIN_DECLS

#ifdef _MSC_VER
#include <math.h>

#define snprintf _snprintf
#define inline  __inline

#ifndef isfinite
#  define isfinite(x) _finite(x)
#endif

extern double _plfit_fmin(double a, double b);
extern double _plfit_round(double x);

#define fmin _plfit_fmin
#define round _plfit_round

#endif /* _MSC_VER */

#ifndef isnan
#  define isnan(x) ((x) != (x))
#endif

#ifndef INFINITY
#  define INFINITY (1.0/0.0)
#endif

#ifndef NAN
#  define NAN ((double)0.0 / (double)DBL_MIN)
#endif

__END_DECLS

#endif /* __PLATFORM_H__ */
