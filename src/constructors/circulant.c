/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_constructors.h"
#include "igraph_interface.h"

// TODO: add documentation

igraph_error_t igraph_circulant(igraph_t *graph, igraph_integer_t n, const igraph_vector_int_t *l, igraph_bool_t directed) {

    igraph_vector_int_t edges;
    igraph_vector_bool_t offset_seen;
    igraph_integer_t i, j;

    if (n < 0) {
        IGRAPH_ERRORF("n = %" IGRAPH_PRId " must be at least 1.", IGRAPH_EINVAL, n);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&offset_seen, n);

    for (i = 0; i < igraph_vector_int_size(l); i++) {
        /* simplify the offset */
        igraph_integer_t offset = VECTOR(*l)[i] % n;
        if (offset < 0) {
            offset += n;
        }
        if (!directed) {
            if (offset >= (n + 1) / 2) {
                offset = n - offset;
            }
        }

        /* only use offset if non-zero and we haven't seen it before */
        if (offset != 0 && !VECTOR(offset_seen)[offset]) {
            for (j = 0; j < n; j++) {
                igraph_vector_int_push_back(&edges, j);
                igraph_vector_int_push_back(&edges, (j + offset) % n);
            }

            VECTOR(offset_seen)[offset] = 1;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    igraph_vector_bool_destroy(&offset_seen);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
