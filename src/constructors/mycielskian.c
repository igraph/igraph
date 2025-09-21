/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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
#include "igraph_operators.h"

#include "igraph_conversion.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

/**
 * \function igraph_mycielskian
 * \brief Generate the Mycielskian of a graph with \p k iterations.
 *
 * \experimental
 *
 * The Mycielskian of a graph is a larger graph formed using a construction due
 * to Jan Mycielski that increases the chromatic number by one while preserving
 * the triangle-free property. The Mycielski construction can be used to create
 * triangle-free graphs with an arbitrarily large chromatic number.
 *
 * </para><para>
 * Let the \c n vertices of the given graph \c G be \c v_1, ..., \c v_n.
 * The Mycielskian of \c G, denoted <code>M(G)</code>, contains \c G itself as
 * a subgraph, together with \c n+1 additional vertices:
 *
 * \ilist
 * \ili A vertex \c u_i corresponding to each vertex \c v_i of \c G.
 * \ili An extra vertex \c w.
 * \endilist
 *
 * </para><para>
 * The edges are added as follows:
 *
 * \ilist
 * \ili Each vertex \c u_i is connected to \c w, forming a star.
 * \ili For each edge <code>(v_i, v_j)</code> in \c G, two new edges are added:
 *   <code>(u_i, v_j)</code> and <code>(v_i, u_j)</code>.
 * \endilist
 *
 * Thus, if \c G has \c n vertices and \c m edges, the Mycielskian <code>M(G)</code>
 * has <code>2n + 1</code> vertices, and <code>3m + n</code> edges.
 *
 * </para><para>
 * igraph uses an alternative construction in two special cases:
 *
 * \ilist
 * \ili The Mycielskian of the null graph is the singleton graph.
 * \ili The Mycielskian of the singleton graph is the two-path.
 * \endilist
 *
 * This ensures that iterative applications of the construction, starting from
 * the null or singleton graph, always yields connected graphs. In fact these
 * are the Mycielski graphs that \ref igraph_mycielski_graph() produces.
 *
 * </para><para>
 * igraph extends the construction to directed graphs, as well as to non-simple
 * graphs, by following the above constructions rules literally.
 *
 * </para><para>
 * This function can apply the Mycielski transformation an arbitrary number of
 * times, controlled by the parameter \p k. The k-th iterated Mycielskian has
 * <code>n_k = (n + 1) * 2^k - 1</code>
 * vertices and
 * <code>m_k = ((2m + 2n + 1) * 3^k - n_{k+1}) / 2</code>
 * edges, where \c n and \c m are the vertex and edge count of the original
 * graph, respectively.
 *
 * \param graph Pointer to the input graph.
 * \param res Pointer to an uninitialized graph object where the Mycielskian
 *        of the input graph will be stored.
 * \param k Integer, the number of Mycielskian iterations (must be non-negative).
 * \return Error code.
 *
 * \sa \ref igraph_mycielski_graph().
 *
 * Time complexity: O(|V| 2^k + |E| 3^k) where |V| and |E| are the vertex and
 * edge counts, respectively.
 */
igraph_error_t igraph_mycielskian(const igraph_t *graph, igraph_t *res, igraph_int_t k) {
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_vector_int_t edges;

    if (k < 0) {
        IGRAPH_ERROR("The number of Mycielski iterations must not be negative.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, /* bycol */ false));

    /* Special case: null graph.
     * We add a single vertex. */
    if (vcount == 0 && k > 0) {
        vcount += 1;
        k--;
    }

    /* Special case: singleton graph.
     * We add a new vertex and connect it to the existing vertex. */
    if (vcount == 1 && k > 0) {
        vcount += 1;
        ecount += 1;
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, 0));
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, 1));
        k--;
    }

    /* No more special cases, we are ready for performing the remaining Mycielski iterations. */

    /* Compute the number of vertices and edges. Since these are exponential in k,
     * overflow checks are important. */

    igraph_int_t new_vcount = vcount;
    igraph_int_t new_ecount = ecount;

    for (igraph_int_t i = 0; i < k; i++) {
        // new edges = 3 * old edges + old vertices
        IGRAPH_SAFE_MULT(new_ecount, 3, &new_ecount);
        IGRAPH_SAFE_ADD(new_ecount, new_vcount, &new_ecount);

        // new vertices = 2 * old vertices + 1
        IGRAPH_SAFE_MULT(new_vcount, 2, &new_vcount);
        IGRAPH_SAFE_ADD(new_vcount, 1, &new_vcount);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(&edges, new_ecount * 2));

    igraph_int_t edge_index = 2 * ecount;  // Current last edge index in edge vector
    igraph_int_t offset = vcount;          // Tracks where new vertices start

    for (igraph_int_t i = 0; i < k; i++) {
        igraph_int_t prev_vcount = offset;  // Number of vertices before this step
        igraph_int_t w = offset * 2;        // The new 'w' node index
        igraph_int_t last_edge_index = edge_index;  // Mark where edges before this step end

        // For each edge before this step, add two new edges
        for (igraph_int_t j = 0; j < last_edge_index; j += 2) {
            igraph_int_t v1 = VECTOR(edges)[j];
            igraph_int_t v2 = VECTOR(edges)[j + 1];

            VECTOR(edges)[edge_index++] = v1;
            VECTOR(edges)[edge_index++] = offset + v2;

            VECTOR(edges)[edge_index++] = v2;
            VECTOR(edges)[edge_index++] = offset + v1;
        }

        // Add edges connecting each `ui` to `w` (forming a star)
        for (igraph_int_t j = prev_vcount; j < w; j++) {
            VECTOR(edges)[edge_index++] = j;
            VECTOR(edges)[edge_index++] = w;
        }

        // Update offset for next step
        offset = offset * 2 + 1;

        IGRAPH_ALLOW_INTERRUPTION();
    }

    IGRAPH_CHECK(igraph_create(res, &edges, new_vcount, igraph_is_directed(graph)));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_mycielski_graph
 * \brief The Mycielski graph of order \p k.
 *
 * The Mycielski graph of order \p k, denoted \c M_k, is a triangle-free graph on
 * \p k vertices with chromatic number \p k. It is defined through the Mycielski
 * construction described in the documentation of \ref igraph_mycielskian().
 *
 * </para><para>
 * Some authors define Mycielski graphs only for <code>k > 1</code>.
 * igraph extends this to all <code>k >= 0</code>.
 * The first few Mycielski graphs are:
 * \olist
 * \oli M_0: Null graph
 * \oli M_1: Single vertex
 * \oli M_2: Path graph with 2 vertices
 * \oli M_3: Cycle graph with 5 vertices
 * \oli M_4: Grötzsch graph (a triangle-free graph with chromatic number 4)
 * \endolist
 *
 * The vertex count of \c M_k is
 * <code>n_k = 3 * 2^(k-2) - 1</code> for <code>k > 1</code> and \c k otherwise.
 * The edge count is
 * <code>m_k = (7 * 3^(k-2) + 1) / 2 - 3 * 2^(k - 2)</code> for <code>k > 1</code>
 * and 0 otherwise.
 *
 * \param graph Pointer to an uninitialized graph object. The generated
 *        Mycielski graph will be stored here.
 * \param k Integer, the order of the Mycielski graph (must be non-negative).
 * \return Error code.
 *
 * \sa \ref igraph_mycielskian().
 *
 * Time complexity: O(3^k), i.e. exponential in \p k.
 */
igraph_error_t igraph_mycielski_graph(igraph_t *graph, igraph_int_t k) {
    igraph_t g;

    if (k < 0) {
        IGRAPH_ERROR("The Mycielski graph order must not be negative.", IGRAPH_EINVAL);
    }

    /* Special cases: M_0 and M_1 are the null graph and singleton graph, respectively. */
    if (k <= 1) {
        IGRAPH_CHECK(igraph_empty(graph, k, IGRAPH_UNDIRECTED));
        return IGRAPH_SUCCESS;
    }

    /* For k >= 2, start construction from P_2. */
    IGRAPH_CHECK(igraph_ring(&g, 2, IGRAPH_UNDIRECTED, false, false));
    IGRAPH_FINALLY(igraph_destroy, &g);

    IGRAPH_CHECK(igraph_mycielskian(&g, graph, k - 2));

    igraph_destroy(&g);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
