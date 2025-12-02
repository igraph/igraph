/*
   igraph library.
   Copyright (C) 2005-2025 The igraph development team

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

#include "igraph_constructors.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"
#include <math.h>

/**
 * \function igraph_hamming_graph
 * \brief Generate a Hamming graph H(d, q)
 *
 * A Hamming graph H(d, q) has \c q^d vertices corresponding to all
 * strings of length \c d over an alphabet of size \c q.
 * Two vertices are adjacent if they differ in exactly one position.
 *
 * \param graph Pointer to an uninitialized graph object; the result will be stored here.
 * \param d Dimension (must be positive)
 * \param q Alphabet size (must be greater than 1)
 * \return Error code.
 *
 * Time complexity: O(q^{2d} * d) in worst case, suitable only for small d.
 */
igraph_error_t igraph_hamming_graph(igraph_t *graph, igraph_int_t d, igraph_int_t q) {

    if (d <= 0 || q <= 1) {
        IGRAPH_ERROR("Both d and q must be positive and q > 1", IGRAPH_EINVAL);
    }

    /* Number of vertices: q^d */
    igraph_real_t n_real = pow(q, d);
    igraph_int_t n = n_real;
    if (n != n_real) {
        IGRAPH_ERRORF("Parameters (%" IGRAPH_PRId ", %" IGRAPH_PRId ") too large for Hamming graph.",
                      IGRAPH_EINVAL, d, q);
    }

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    /* Represent vertices as numbers 0...(q^d - 1) */
    for (igraph_int_t i = 0; i < n; i++) {
        for (igraph_int_t j = i + 1; j < n; j++) {
            igraph_int_t diff = 0;
            igraph_int_t a = i, b = j;

            /* Compare base-q digits */
            for (igraph_int_t k = 0; k < d; k++) {
                if ((a % q) != (b % q)) diff++;
                a /= q;
                b /= q;
            }

            if (diff == 1) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            }
        }
        IGRAPH_ALLOW_INTERRUPTION();
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
