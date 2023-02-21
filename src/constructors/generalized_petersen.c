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

#include "math/safe_intop.h"

/**
 * \function igraph_generalized_petersen
 * \brief Creates a Generalized Petersen graph.
 *
 * The generalized Petersen graph <code>G(n, k)</code> consists of \p n vertices
 * \c v_0, ..., \c v_n forming an "outer" cycle graph, and \p n additional vertices
 * \c u_0, ..., \c u_n forming an "inner" circulant graph where <code>u_i</code>
 * is connected to <code>u_(i + k mod n)</code>. Additionally, all \c v_i are
 * connected to \c u_i.
 *
 * </para><para>
 * <code>G(n, k)</code> has \c 2n vertices and \c 3n edges. The Petersen graph
 * itself is <code>G(5, 2)</code>.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * M. E. Watkins,
 * A Theorem on Tait Colorings with an Application to the Generalized Petersen Graphs,
 * Journal of Combinatorial Theory 6, 152-164 (1969).
 * https://doi.org/10.1016%2FS0021-9800%2869%2980116-X
 *
 * \param graph Pointer to an uninitialized graph object, the result will
 * be stored here.
 * \param n Integer, \c n is the number of vertices in the inner and outer
 * cycle/circulant graphs. It must be at least 3.
 * \param k Integer, \c k is the shift of the circulant graph. It must be
 * positive and less than <code>n/2</code>.
 * \return Error code.
 *
 * \sa \ref igraph_famous() for the original Petersen graph.
 *
 * Time complexity: O(|V|), the number of vertices in the graph.
 */
igraph_error_t igraph_generalized_petersen(igraph_t *graph, igraph_integer_t n, igraph_integer_t k) {
    /* This is a generalized Petersen graph constructor */
    igraph_vector_int_t edges;
    igraph_integer_t no_of_nodes, no_of_edges2;
    igraph_integer_t i;

    if (n < 3) {
        IGRAPH_ERRORF("n = %" IGRAPH_PRId " must be at least 3.", IGRAPH_EINVAL, n);
    }

    IGRAPH_SAFE_MULT(n, 2, &no_of_nodes);

    /* The seemingly redundant k < n check avoids integer overflow on 2*k in 2*k < n.
     * Note that 2*n has already been checked not to overflow above. */
    if (! (k > 0 && k < n && 2*k < n)) {
        IGRAPH_ERRORF("k = %" IGRAPH_PRId " must be positive and less than n/2 with n = %" IGRAPH_PRId ".", IGRAPH_EINVAL, k, n);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_SAFE_MULT(n, 6, &no_of_edges2);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));

    for (i = 0; i < n; i++) {
        igraph_vector_int_push_back(&edges, i);
        igraph_vector_int_push_back(&edges, (i + 1) % n);
        igraph_vector_int_push_back(&edges, i);
        igraph_vector_int_push_back(&edges, i + n);
        igraph_vector_int_push_back(&edges, i + n);
        igraph_vector_int_push_back(&edges, ((i + k) % n) + n);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
