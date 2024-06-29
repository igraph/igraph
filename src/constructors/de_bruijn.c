/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

/**
 * \function igraph_de_bruijn
 * \brief Generate a de Bruijn graph.
 *
 * A de Bruijn graph represents relationships between strings. An alphabet
 * of \c m letters are used and strings of length \c n are considered.
 * A vertex corresponds to every possible string and there is a directed edge
 * from vertex \c v to vertex \c w if the string of \c v can be transformed into
 * the string of \c w by removing its first letter and appending a letter to it.
 *
 * </para><para>
 * Please note that the graph will have \c m to the power \c n vertices and
 * even more edges, so probably you don't want to supply too big numbers for
 * \c m and \c n.
 *
 * </para><para>
 * De Bruijn graphs have some interesting properties, please see another source,
 * e.g. Wikipedia for details.
 *
 * \param graph Pointer to an uninitialized graph object, the result will be
 *        stored here.
 * \param m Integer, the number of letters in the alphabet.
 * \param n Integer, the length of the strings.
 * \return Error code.
 *
 * \sa \ref igraph_kautz().
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of edges.
 */
igraph_error_t igraph_de_bruijn(igraph_t *graph, igraph_integer_t m, igraph_integer_t n) {

    /* m - number of symbols */
    /* n - length of strings */

    igraph_integer_t no_of_nodes, no_of_edges;
    igraph_vector_int_t edges;
    igraph_integer_t i, j;
    int iter = 0;

    if (m < 0 || n < 0) {
        IGRAPH_ERROR("`m' and `n' should be non-negative in a de Bruijn graph",
                     IGRAPH_EINVAL);
    }

    if (n == 0) {
        return igraph_empty(graph, 1, IGRAPH_DIRECTED);
    }
    if (m == 0) {
        return igraph_empty(graph, 0, IGRAPH_DIRECTED);
    }

    {
        igraph_real_t no_of_nodes_real = pow(m, n);
        no_of_nodes = no_of_nodes_real;
        if (no_of_nodes != no_of_nodes_real) {
            IGRAPH_ERRORF("Parameters (%" IGRAPH_PRId ", %" IGRAPH_PRId ") too large for De Bruijn graph.", IGRAPH_EINVAL,
                          m, n);
        }
    }
    /* no_of_edges = m * no_of_nodes */
    IGRAPH_SAFE_MULT(no_of_nodes, m, &no_of_edges);

    {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges2);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
    }

    for (i = 0; i < no_of_nodes; i++) {
        igraph_integer_t basis = (i * m) % no_of_nodes;
        for (j = 0; j < m; j++) {
            igraph_vector_int_push_back(&edges, i);
            igraph_vector_int_push_back(&edges, basis + j);
        }
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_DIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
