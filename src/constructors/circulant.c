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

/**
 * \function igraph_circulant
 * \brief Creates a circulant graph.
 *
 * A circulant graph <code>G(n, l)</code> consists of \p n vertices \c v_0, ...,
 * \c v_(n-1) such that for each \c l_i in the list of offsets \p l, \c v_j is
 * connected to <code> v_((j + l_i) mod n) </code> for all j.
 *
 * </para><para>
 * The function works with both directed and undirected graphs. Multiple edges are
 * merged, and self loops are ignored. All offsets are taken modulo \p n, so values greater than n as well as negative offsets are permitted.
 *
 * \param graph Pointer to an uninitialized graph object, the result will
 * be stored here.
 * \param n Integer, \p n is the number of vertices in the circulant graph. It must
 * be at least 1.
 * \param l Integer vector, \p l is a list of the offsets within the circulant graph.
 * \param directed Boolean, \p directed determines whether the graph should be directed.
 * \return Error code.
 *
 * \sa \ref igraph_ring(), \ref igraph_generalized_petersen().
 *
 * Time complexity: O(|V||L|), the number of vertices in the graph times the number
 * of offsets.
 */

igraph_error_t igraph_circulant(igraph_t *graph, igraph_integer_t n, const igraph_vector_int_t *l, igraph_bool_t directed) {

    igraph_vector_int_t edges;
    igraph_vector_bool_t offset_seen;
    igraph_integer_t i, j;
    igraph_integer_t limit;

    if (n <= 0) {
        IGRAPH_ERRORF("number of nodes = %" IGRAPH_PRId " must be at least 1.", IGRAPH_EINVAL, n);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&offset_seen, n);

    VECTOR(offset_seen)[0] = 1; /* do not allow self loops */

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
        if (!VECTOR(offset_seen)[offset]) {
            if (n % 2 == 0 && offset == n / 2 && !directed) {
                limit = n / 2; /* this to avoid doubling up the n/2 offset for even n and undirected graph */
            } else {
                limit = n;
            }
            for (j = 0; j < limit; j++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, (j + offset) % n));
            }

            VECTOR(offset_seen)[offset] = 1;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    igraph_vector_bool_destroy(&offset_seen);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
