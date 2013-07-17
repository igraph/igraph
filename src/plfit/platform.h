/* platform.h
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
#define snprintf sprintf_s
#define inline  __inline
#define isnan(x) _isnan(x)
#define isfinite(x) _finite(x)
#endif

#ifndef INFINITY
#  define INFINITY (1.0/0.0)
#endif

#ifndef NAN
#  define NAN (INFINITY-INFINITY)
#endif

__END_DECLS

#endif /* __PLATFORM_H__ */
