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

#ifndef IGRAPH_HACKS_INTERNAL_H
#define IGRAPH_HACKS_INTERNAL_H

#include "config.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef HAVE_STRDUP
    #define strdup igraph_i_strdup
    char* igraph_i_strdup(const char *s);
#endif

#ifndef HAVE_STPCPY
    #define stpcpy igraph_i_stpcpy
    char* igraph_i_stpcpy(char* s1, const char* s2);
#else
    #ifndef HAVE_STPCPY_SIGNATURE
        char* stpcpy(char* s1, const char* s2);
    #endif
#endif

__END_DECLS

#endif
