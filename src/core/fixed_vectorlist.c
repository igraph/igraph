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

#include "core/fixed_vectorlist.h"

void igraph_fixed_vectorlist_destroy(igraph_fixed_vectorlist_t *l) {
    igraph_vector_int_list_destroy(&l->vecs);
}

igraph_error_t igraph_fixed_vectorlist_convert(
    igraph_fixed_vectorlist_t *l, const igraph_vector_int_t *from,
    igraph_integer_t size
) {
    igraph_vector_int_t sizes;
    igraph_integer_t i, no = igraph_vector_int_size(from), to;

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&l->vecs, size);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, size);

    for (i = 0; i < no; i++) {
        to = VECTOR(*from)[i];
        if (to >= 0) {
            VECTOR(sizes)[to] += 1;
        }
    }

    for (i = 0; i < no; i++) {
        to = VECTOR(*from)[i];
        if (to >= 0) {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&l->vecs, to);
            IGRAPH_CHECK(igraph_vector_int_push_back(v, i));
        }
    }

    igraph_vector_int_destroy(&sizes);
    IGRAPH_FINALLY_CLEAN(2);  /* + l->vecs */

    return IGRAPH_SUCCESS;
}
