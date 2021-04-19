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

#ifndef IGRAPH_MEMORY_H
#define IGRAPH_MEMORY_H

#include <stdlib.h>
#include "igraph_decls.h"

__BEGIN_DECLS

#define IGRAPH_CALLOC(n,t)    (t*) calloc( (n) > 0 ? (size_t)(n) : (size_t)1, sizeof(t) )
#define IGRAPH_REALLOC(p,n,t) (t*) realloc((void*)(p), (n) > 0 ? (size_t)((n)*sizeof(t)) : (size_t)1)
#define IGRAPH_FREE(p)        (free( (void *)(p) ), (p) = NULL)

#define igraph_Calloc IGRAPH_CALLOC
#define igraph_Realloc IGRAPH_REALLOC
#define igraph_Free IGRAPH_FREE

IGRAPH_EXPORT void igraph_free(void *p);
IGRAPH_EXPORT void *igraph_malloc(size_t n);

__END_DECLS

#endif
