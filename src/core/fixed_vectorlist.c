/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_memory.h"

#include "core/fixed_vectorlist.h"

void igraph_fixed_vectorlist_destroy(igraph_fixed_vectorlist_t *l) {
    long int i, n = igraph_vector_ptr_size(&l->v);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(l->v)[i];
        if (v) {
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_destroy(&l->v);
    igraph_free(l->vecs);
}

int igraph_fixed_vectorlist_convert(igraph_fixed_vectorlist_t *l,
                                    const igraph_vector_t *from,
                                    long int size) {

    igraph_vector_t sizes;
    long int i, no = igraph_vector_size(from);

    l->vecs = IGRAPH_CALLOC(size, igraph_vector_t);
    if (!l->vecs) {
        IGRAPH_ERROR("Cannot merge attributes for simplify",
                     IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, l->vecs);
    IGRAPH_CHECK(igraph_vector_ptr_init(&l->v, size));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &l->v);
    IGRAPH_VECTOR_INIT_FINALLY(&sizes, size);

    for (i = 0; i < no; i++) {
        long int to = (long int) VECTOR(*from)[i];
        if (to >= 0) {
            VECTOR(sizes)[to] += 1;
        }
    }
    for (i = 0; i < size; i++) {
        igraph_vector_t *v = &(l->vecs[i]);
        IGRAPH_CHECK(igraph_vector_init(v, (long int) VECTOR(sizes)[i]));
        igraph_vector_clear(v);
        VECTOR(l->v)[i] = v;
    }
    for (i = 0; i < no; i++) {
        long int to = (long int) VECTOR(*from)[i];
        if (to >= 0) {
            igraph_vector_t *v = &(l->vecs[to]);
            igraph_vector_push_back(v, i);
        }
    }

    igraph_vector_destroy(&sizes);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}
