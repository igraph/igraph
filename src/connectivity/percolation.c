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

#include "igraph_components.h"

#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_types.h"
#include "igraph_vector.h"

/**
 * \function igraph_i_percolate_edge
 * \brief Percolates a single edge.
 * \param links Vector representing parents.
 * \param sizes sizes[i] is the number of children of links[i]
 * \param biggest The biggest value in sizes, is updated if a bigger cluster is created.
 * \param a A vertex incident to the edge.
 * \param b The other vertex incident to the edge.
 */

static void percolate_edge(igraph_vector_int_t *links,
                           igraph_vector_int_t *sizes,
                           igraph_integer_t *biggest,
                           igraph_integer_t a,
                           igraph_integer_t b) {

    // Find head of each tree
    while (VECTOR(*links)[a] != a) {
        VECTOR(*links)[a] = VECTOR(*links)[VECTOR(*links)[a]];
        a = VECTOR(*links)[a];
    }
    while (VECTOR(*links)[b] != b) {
        VECTOR(*links)[b] = VECTOR(*links)[VECTOR(*links)[b]];
        b = VECTOR(*links)[b];
    }

    // If they are already connected, exit early
    if (a == b) {
        return;
    }

    // Make smaller child of larger
    igraph_integer_t parent, child;
    if (VECTOR(*sizes)[a] < VECTOR(*sizes)[b]) {
        parent = b;
        child = a;
    } else {
        parent = a;
        child = b;
    }

    // Make a child of b
    VECTOR(*links)[child] = parent;
    VECTOR(*sizes)[parent] += VECTOR(*sizes)[child];

    // If made new biggest component, update biggest
    if (VECTOR(*sizes)[parent] >= *biggest) {
        *biggest = VECTOR(*sizes)[parent];
    }
}

/**
 * \function igraph_edgelist_percolation
 * \brief Gives the size of the largest connected component as edges are added.
 *
 * \experimental
 *
 * Calculates the size of the largest component as edges are added to a graph in the order given.
 * If the edge-sequence and output are reversed, they are the size of the largest component as
 * edges are removed.
 *
 * \param edges Vector of edges, where the n-th edge has endpoints edges[2n] and edges [2n+1].
 * \param output output[n] is the size of the largest connected component after edge n is added.
 *        Will be resized.
 *
 * \return Error code
 *
 * Time complexity: O(|E| * a(|E|)) where a is the inverse Ackermann function,
 *                  for all practical purposes it is not above 5.
 */

igraph_error_t igraph_edgelist_percolation(
        const igraph_vector_int_t *edges,
        igraph_vector_int_t* output) {

    igraph_integer_t biggest = 0;
    igraph_integer_t lower, upper;

    // Handle edge case of no edges.
    if (igraph_vector_int_size(edges) == 0) {
        IGRAPH_CHECK(igraph_vector_int_resize(output, 0));
        return IGRAPH_SUCCESS;
    }
    igraph_vector_int_minmax(edges, &lower, &upper);

    if (lower < 0) {
        IGRAPH_ERROR("Invalid (negative) vertex index.", IGRAPH_EINVVID);
    }

    igraph_vector_int_t sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, upper + 1);

    igraph_vector_int_t links;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&links, upper + 1);

    for (igraph_integer_t i = 0; i < upper + 1; i++) {
        VECTOR(sizes)[i] = 1;
        VECTOR(links)[i] = i;
    }

    igraph_integer_t edge_count = igraph_vector_int_size(edges);
    if (edge_count % 2 == 1) {
        IGRAPH_ERROR("Invalid edge list, odd number of elements.", IGRAPH_EINVAL);
    }
    edge_count >>= 1;
    IGRAPH_CHECK(igraph_vector_int_resize(output, edge_count));

    for (igraph_integer_t i = 0; i < edge_count; i++) {
        percolate_edge(&links, &sizes, &biggest, VECTOR(*edges)[2 * i], VECTOR(*edges)[2 * i + 1]);
        VECTOR(*output)[i] = biggest;
    }

    igraph_vector_int_destroy(&sizes);
    igraph_vector_int_destroy(&links);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bond_percolation
 * \brief Calculates the bond-percolation curve of a graph.
 *
 * \experimental
 *
 * Calculates the bond-percolation curve, or the size of the largest component as
 * edges are added to the graph in the order given. If both the output and the input
 * are reversed, it is the size of the largest component as edges are removed.
 *
 * \param graph The input graph
 * \param output output[i] is the size of the largest component after adding
 *    <code>i+1</code> edges, will be created.
 * \param edge_order The order of the edges, will be generated at random if \c NULL.
 * \return Error code.
 *
 * Time complexity: O(|V| + |E| * a(|E|)) where a is the inverse Ackermann function,
 *                  for all practical purposes it is not above 5.
 */

igraph_error_t igraph_bond_percolation(
        const igraph_t *graph,
        igraph_vector_int_t *output,
        const igraph_vector_int_t *edge_order) {

    const igraph_vector_int_t *p_edge_order;
    igraph_vector_int_t i_edge_order;
    igraph_vector_int_t edges;

    if (edge_order == NULL) {
        // Use random edge order

        IGRAPH_CHECK(igraph_vector_int_init_range(&i_edge_order, 0, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &i_edge_order);
        IGRAPH_CHECK(igraph_vector_int_shuffle(&i_edge_order));

        p_edge_order = &i_edge_order;
    } else {
        p_edge_order = edge_order;
    }

    IGRAPH_CHECK(igraph_vector_int_init(&edges, 2 * igraph_vector_int_size(p_edge_order)));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &edges);

    IGRAPH_CHECK(igraph_edges(graph, igraph_ess_vector(p_edge_order), &edges));
    IGRAPH_CHECK(igraph_edgelist_percolation(&edges, output));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    if (edge_order == NULL) {
        igraph_vector_int_destroy(&i_edge_order);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t percolate_site(const igraph_t *graph,
                                     igraph_vector_int_t *links,
                                     igraph_vector_int_t *sizes,
                                     igraph_integer_t *biggest,
                                     igraph_integer_t vertex,
                                     igraph_vector_int_t *neighbors) {

    if (VECTOR(*sizes)[vertex] != 0) {
        IGRAPH_ERROR("Duplicate vertex in vertex order vector.", IGRAPH_EINVAL);
    }

    VECTOR(*sizes)[vertex] = 1;

    igraph_integer_t neighbor_count;

    IGRAPH_CHECK(igraph_neighbors(graph, neighbors, vertex, IGRAPH_IN));

    neighbor_count = igraph_vector_int_size(neighbors);
    for (igraph_integer_t i = 0; i < neighbor_count; i++) {
        if (VECTOR(*sizes)[VECTOR(*neighbors)[i]] == 0) {
            continue;
        }
        percolate_edge(links, sizes, biggest, vertex, VECTOR(*neighbors)[i]);
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_site_percolation
 * \brief The size of the largest component as vertices are added to a graph.
 *
 * \experimental
 *
 * Calculates the percolation curve, or the size of the largest connected component as
 * vertices are added in the order given. If the output is reversed, it is the size of
 * the largest component as vertices are removed in the reverse of the order given.
 * If there is no vertex order given, it will generate one.
 *
 * \param graph The graph that vertices are assumed to be in.
 * \param output <code>output[i]</code> will contain the size of the largest component
 *        after adding <code>vertex_order[i]</code>. Will be resized.
 * \param vertex_order The order the vertices will be added in.
 *        Will raise error if there are duplicates, or a vertex is missing.
 *        If \c NULL, a random order will be used.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_EINVAL
 *             The vertex list is invalid, see above.
 *        \endclist
 *
 * Time Complexity: O(|V| + |E| * a(|E|)) where a is the inverse Ackermann function,
 *                  for all practical purposes it is not above 5.
 */

igraph_error_t igraph_site_percolation(
        const igraph_t *graph,
        igraph_vector_int_t *output,
        const igraph_vector_int_t *vertex_order) {

    const igraph_vector_int_t *p_vertex_order;
    igraph_vector_int_t i_vertex_order;
    const igraph_integer_t vcount = igraph_vcount(graph);

    if (vertex_order == NULL) {
        // Use random vertex order

        IGRAPH_CHECK(igraph_vector_int_init_range(&i_vertex_order, 0, vcount));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &i_vertex_order);
        IGRAPH_CHECK(igraph_vector_int_shuffle(&i_vertex_order));

        p_vertex_order = &i_vertex_order;
    } else {
        p_vertex_order = vertex_order;
    }

    // Initialize variables
    if (igraph_vector_int_size(p_vertex_order) != vcount) {
        IGRAPH_ERRORF("Vertex order vector length (%" IGRAPH_PRId") "
                      "does not match the vertex count (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_int_size(p_vertex_order), vcount);
    }
    igraph_integer_t biggest = 1;

    igraph_vector_int_t sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, vcount);

    igraph_vector_int_t links;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&links, vcount);

    IGRAPH_CHECK(igraph_vector_int_resize(output, vcount));

    for (igraph_integer_t i = 0; i < vcount; i++) {
        VECTOR(sizes)[i] = 0;
        VECTOR(links)[i] = i;
    }

    igraph_vector_int_t neighbors;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 0);
    for (igraph_integer_t i = 0; i < vcount; i++) {
        IGRAPH_CHECK(percolate_site(graph, &links, &sizes, &biggest, VECTOR(*p_vertex_order)[i], &neighbors));
        VECTOR(*output)[i] = biggest;
    }
    igraph_vector_int_destroy(&neighbors);
    IGRAPH_FINALLY_CLEAN(1);

    for (igraph_integer_t i = 0; i < vcount; i++) {
        if (VECTOR(sizes)[i] == 0) {
            IGRAPH_ERROR("Vertex order vector is missing vertices from graph.", IGRAPH_EINVAL);
        }
    }

    // Cleanup
    igraph_vector_int_destroy(&links);
    igraph_vector_int_destroy(&sizes);
    IGRAPH_FINALLY_CLEAN(2);

    if (vertex_order == NULL) {
        igraph_vector_int_destroy(&i_vertex_order);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
