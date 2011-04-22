/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_THREADING_H
#define IGRAPH_THREADING_H

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

/* TODO: we are detecting support for thread-local storage here based on
 * the compiler and not the actual capabilities of the system; for instance,
 * it could happen that the compiler is gcc but it does not support
 * thread-local storage on the current platform. We should update the
 * ./configure script to detect TLS-support properly if this causes problems.
 */

#ifdef HAVE_TLS

#if defined(__GNUC__) || defined(__SUNPRO_C) || defined(__BORLANDC__)
   /* gcc, Sun Studio, Intel C++ on Linux, Borland C */
#  define IGRAPH_THREAD_LOCAL __thread
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
   /* Microsoft Visual C++, Intel C++ on Windows */
#  define IGRAPH_THREAD_LOCAL __declspec(thread)
#else
   /* Other compilers -- play it safe */
#  define IGRAPH_THREAD_LOCAL
#  warning "Thread-local storage is not supported by this compiler."
#endif

#else

#define IGRAPH_THREAD_LOCAL

#endif  /* HAVE_TLS */

__END_DECLS

#endif

