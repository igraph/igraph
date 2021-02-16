/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_TYPES_INTERNAL_H
#define IGRAPH_TYPES_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Vectorlist, fixed length                           */
/* -------------------------------------------------- */

typedef struct igraph_fixed_vectorlist_t {
    igraph_vector_t *vecs;
    igraph_vector_ptr_t v;
    long int length;
} igraph_fixed_vectorlist_t;

void igraph_fixed_vectorlist_destroy(igraph_fixed_vectorlist_t *l);
int igraph_fixed_vectorlist_convert(igraph_fixed_vectorlist_t *l,
                                    const igraph_vector_t *from,
                                    long int size);

__END_DECLS

#endif
