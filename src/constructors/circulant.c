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

#include "math/safe_intop.h"

/**
 * \function igraph_circulant
 * \brief Creates a circulant graph.
 *
 * A circulant graph <code>G(n, shifts)</code> consists of \p n vertices <code>v_0</code>, ...,
 * <code>v_(n-1)</code> such that for each \c s_i in the list of offsets \p shifts, \c v_j is
 * connected to <code>v_((j + s_i) mod n)</code> for all j.
 *
 * </para><para>
 * The function can generate either directed or undirected graphs. It does not generate
 * multi-edges or self-loops.
 *
 * \param graph Pointer to an uninitialized graph object, the result will
 *   be stored here.
 * \param n Integer, the number of vertices in the circulant graph.
 * \param shifts Integer vector, a list of the offsets within the circulant graph.
 * \param directed Boolean, whether to create a directed graph.
 * \return Error code.
 *
 * \sa \ref igraph_ring(), \ref igraph_generalized_petersen(), \ref igraph_extended_chordal_ring()
 *
 * Time complexity: O(|V||shifts|), the number of vertices in the graph times the number
 * of shifts.
 */

igraph_error_t igraph_circulant(igraph_t *graph, igraph_integer_t n, const igraph_vector_int_t *shifts, igraph_bool_t directed) {

    igraph_vector_int_t edges;
    igraph_vector_bool_t shift_seen;
    igraph_integer_t i, j;
    igraph_integer_t limit;
    igraph_integer_t shift_size = igraph_vector_int_size(shifts);

    if (n < 0) {
        IGRAPH_ERRORF("Number of nodes = %" IGRAPH_PRId " must be non-negative.", IGRAPH_EINVAL, n);
    }
    if (n == 0) {
        return igraph_empty(graph, 0, directed);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    {
        igraph_integer_t size;
        IGRAPH_SAFE_MULT(n, shift_size, &size);
        IGRAPH_SAFE_MULT(size, 2, &size);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, size));
    }

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&shift_seen, n);
    VECTOR(shift_seen)[0] = true; /* do not allow self loops */

    for (i = 0; i < shift_size; i++) {
        /* simplify the shift */
        igraph_integer_t shift = VECTOR(*shifts)[i] % n;
        if (shift < 0) {
            shift += n;
        }
        if (!directed) {
            if (shift >= (n + 1) / 2) {
                shift = n - shift;
            }
        }

        /* only use shift if non-zero and we haven't seen it before */
        if (!VECTOR(shift_seen)[shift]) {
            if (n % 2 == 0 && shift == n / 2 && !directed) {
                limit = n / 2; /* this to avoid doubling up the n/2 shift for even n and undirected graph */
            } else {
                limit = n;
            }
            for (j = 0; j < limit; j++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, (j + shift) % n));
            }

            VECTOR(shift_seen)[shift] = true;
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    igraph_vector_bool_destroy(&shift_seen);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
