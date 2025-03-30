/*
   IGraph library.
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

#include "igraph_structural.h"

#include "igraph_interface.h"

/* Given a vertex order, computes a mapping such that for each index = vertex ID, the
 * corresponding value is its position in the order.
 *
 * Index: ordering  -> Index: vertex ID
 * Value: vertex ID    Value: ordering
 */
static igraph_error_t get_vertex_order_map(const igraph_vector_int_t *vertex_order,
                                           igraph_vector_int_t *map) {
    igraph_integer_t size = igraph_vector_int_size(vertex_order);
    IGRAPH_CHECK(igraph_vector_int_resize(map, size));
    for (igraph_integer_t i = 0; i < size; i++) {
        VECTOR(*map)[VECTOR(*vertex_order)[i]] = i;
    }
    return IGRAPH_SUCCESS;
}

/* Returns the total number of possible edges given the number of vertices in a graph,
 * whether it is directed, and whether loops should be assumed possible.
 *
 * No loops, undirected: n * (n-1) / 2
 * No loops, directed:   n * (n-1)
 * Loops, undirected:    n * (n+1) / 2
 * Loops, directed:      n^2
 */
static igraph_real_t total_possible_edges(igraph_integer_t numVertices,
                                          igraph_bool_t directed,
                                          igraph_bool_t loops) {
    igraph_real_t nv = (igraph_real_t) numVertices;
    if (!loops) {
        return (directed ? nv * (nv - 1) : nv * (nv - 1) / 2);
    } else {
        return (directed ? nv * nv : nv * (nv + 1) / 2);
    }
}

/**
 * \function igraph_rich_club_density_sequence
 * \brief Density sequence of subgraphs formed by sequential vertex removal.
 *
 * This function takes a graph and a vertex ordering as input, sequentially removes the
 * vertices in the given order, and calculates the density of the subgraph after each
 * removal.
 *
 * </para><para>
 * Density is calculated as the ratio of the number of edges (or total edge weight, if
 * weighted) to the number of total possible edges in the graph. The latter is dependent
 * on whether the graph is directed and whether self-loops are assumed to be possible: for
 * undirected graphs without self-loops, this total is given by <code>n(n-1)/2</code>,
 * and for directed graphs by <code>n(n-1)</code>. When self-loops are allowed, these
 * are adjusted to <code>n(n+1)/2</code> for undirected and <code>n^2</code> for directed
 * graphs.
 *
 * </para><para>
 * Vertex order can be sorted by degree so that the resulting density sequence helps
 * reveal how interconnected a graph is across different degree levels, or the presence
 * of a "rich-club" effect.
 *
 * \param graph The graph object to analyze.
 * \param vertex_order Vector giving the order in which vertices are removed.
 * \param directed If false, directed graphs will be treated as undirected.
 *    Ignored in undirected graphs.
 * \param loops Whether self-loops are assumed to be possible.
 * \param weights Vector with weight of edges, if considered (if not, this should be \c NULL).
 * \param res Integer vector containing the result. It must be initialized and
 *    will be resized to be the appropriate size.
 *
 * \return Error code: \c IGRAPH_EINVAL: invalid vertex_order vector and/or weight vector
 * lengths
 *
 * Time complexity: O(|V| + |E|)
 * where |V| is the number of vertices and |E| the number of edges in the graph given.
 *
 * \sa \ref igraph_density(), which uses the same calculation of total possible edges.
 */
igraph_error_t igraph_rich_club_density_sequence(const igraph_t *graph,
                                                 const igraph_vector_int_t *vertex_order,
                                                 igraph_bool_t directed,
                                                 igraph_bool_t loops,
                                                 igraph_vector_t *weights,
                                                 igraph_vector_t *res) {

    const igraph_integer_t vcount = igraph_vcount(graph);
    const igraph_integer_t ecount = igraph_ecount(graph);

    // Error handling: invalid vertex_order and weights sizes
    if (igraph_vector_int_size(vertex_order) != vcount) {
        IGRAPH_ERRORF("Vertex order vector length (%" IGRAPH_PRId ") does not match "
                      "number of vertices (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_int_size(vertex_order), vcount);
    }
    if (weights && igraph_vector_size(weights) != ecount) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match "
                      "number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), ecount);
    }

    if (! igraph_is_directed(graph)) {
        directed = false;
    }

    igraph_vector_t edges_remaining_after;
    igraph_vector_int_t order_of;

    IGRAPH_VECTOR_INIT_FINALLY(&edges_remaining_after, vcount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order_of, vcount);

    IGRAPH_CHECK(igraph_vector_resize(res, vcount)); // resize res
    IGRAPH_CHECK(get_vertex_order_map(vertex_order, &order_of)); // get map of vertex order

    igraph_bool_t warning_issued = false;

    // edgesRemainingAfter vector: number of edges (or total edge weight) removed by index
    for (igraph_integer_t eid = 0; eid < ecount; eid++) {
        igraph_integer_t v1 = IGRAPH_FROM(graph, eid);
        igraph_integer_t v2 = IGRAPH_TO(graph, eid);
        if (!loops && !warning_issued && v1 == v2) {
            IGRAPH_WARNING("Self-loop encountered while `loops = false`. "
                           "Proceeding as if loops were not possible (density computed accordingly)");
            warning_issued = true;
        }
        igraph_integer_t order_v1 = VECTOR(order_of)[v1]; // order of endpoints
        igraph_integer_t order_v2 = VECTOR(order_of)[v2];

        igraph_integer_t edge_removal_index = (order_v1 < order_v2 ? order_v1 : order_v2);
        VECTOR(edges_remaining_after)[edge_removal_index] += (weights ? VECTOR(*weights)[eid] : 1);
    }

    // edgesRemainingAfter vector: edges (or total edge weight) remaining after i removals
    igraph_real_t total = 0;
    for (igraph_integer_t i = vcount - 1; i >= 0; i--) {
        total += VECTOR(edges_remaining_after)[i];
        VECTOR(edges_remaining_after)[i] = total;
    }

    // density calculation
    for (igraph_integer_t i = 0; i < vcount; i++) {
        // (numVertices - i) = the number of vertices left on this loop
        VECTOR(*res)[i] =
            VECTOR(edges_remaining_after)[i] / total_possible_edges(vcount - i, directed, loops);
    }

    igraph_vector_int_destroy(&order_of);
    igraph_vector_destroy(&edges_remaining_after);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
