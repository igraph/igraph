/*
   igraph library.
   Copyright (C) 2003-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_MEMORY_H
#define IGRAPH_MEMORY_H

#include "igraph_decls.h"

#include <stdint.h>
#include <stdlib.h>

IGRAPH_BEGIN_C_DECLS

/* Helper macro to check if n*sizeof(t) overflows in IGRAPH_CALLOC and IGRAPH_REALLOC */
#define IGRAPH_I_ALLOC_CHECK_OVERFLOW(n,t,expr) \
    (t*) ((0 <= (n) && ((size_t)(n)) <= SIZE_MAX / sizeof(t)) ? (expr) : NULL)

#define IGRAPH_CALLOC(n,t)    IGRAPH_I_ALLOC_CHECK_OVERFLOW(n, t, calloc(sizeof(t) * ((n) > 0 ? (n) : 1), 1))
#define IGRAPH_MALLOC(n)      malloc( (size_t) ((n) > 0 ? (n) : 1) )
#define IGRAPH_REALLOC(p,n,t) IGRAPH_I_ALLOC_CHECK_OVERFLOW(n, t, realloc((void*)(p), sizeof(t) * ((n) > 0 ? (n) : 1)))
#define IGRAPH_FREE(p)        (free( (void *)(p) ), (p) = NULL)

IGRAPH_EXPORT void *igraph_calloc(size_t count, size_t size);
IGRAPH_EXPORT void *igraph_malloc(size_t size);
IGRAPH_EXPORT void *igraph_realloc(void* ptr, size_t size);
IGRAPH_EXPORT void igraph_free(void *ptr);

IGRAPH_END_C_DECLS

#endif
