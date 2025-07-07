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


// test if node is already connected to something
// If it is a root, its size will be greater than 1
// If it is not a root then it will link somewhere.
static bool is_connected(igraph_vector_int_t *links, igraph_vector_int_t *sizes, igraph_integer_t vert) {
    return VECTOR(*links)[vert] != vert || VECTOR(*sizes)[vert] > 1;
}

/**
 * \function igraph_edgelist_percolation
 * \brief The size of the largest connected component as vertex pairs are connected.
 *
 * \experimental
 *
 * Calculates the size of the largest component as edges are added to a graph in
 * the order given. This function differs from \ref igraph_bond_percolation() in that
 * it take a list of vertex pairs as input.
 *
 * \param edges Vector of edges, where the i-th edge has endpoints <code>edges[2i]</code>
 *    and <code>edges[2i+1]</code>.
 * \param giant_size <code>giant_size[i]</code> will contain the size of the largest connected
 *    component after edge \c i is added.
 * \param vertex_count <code>vertex_count[i]</code> will contain the number of vertices
 * with at least one edge after edge \c i is added.
 * \return Error code.
 *
 * \sa \ref igraph_bond_percolation() to specify edges by their ID in a graph object.
 *
 * Time complexity: O(|E| a(|E|)) where a is the inverse Ackermann function,
 * for all practical purposes it is not above 5.
 */

igraph_error_t igraph_edgelist_percolation(
        const igraph_vector_int_t *edges,
        igraph_vector_int_t *giant_size,
        igraph_vector_int_t *vertex_count) {

    igraph_integer_t biggest = 0;
    igraph_integer_t connected = 0;
    igraph_integer_t lower, upper;

    // Handle edge case of no edges.
    if (igraph_vector_int_size(edges) == 0) {
        if (giant_size != NULL) {
            IGRAPH_CHECK(igraph_vector_int_resize(giant_size, 0));
        }
        if (vertex_count != NULL) {
            IGRAPH_CHECK(igraph_vector_int_resize(vertex_count, 0));
        }
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
    IGRAPH_CHECK(igraph_vector_int_resize(giant_size, edge_count));
    IGRAPH_CHECK(igraph_vector_int_resize(vertex_count, edge_count));

    for (igraph_integer_t i = 0; i < edge_count; i++) {
        if (!is_connected(&links, &sizes, VECTOR(*edges)[2 * i]))     {connected++;}
        if (!is_connected(&links, &sizes, VECTOR(*edges)[2 * i+ 1 ])) {connected++;}
        percolate_edge(&links, &sizes, &biggest, VECTOR(*edges)[2 * i], VECTOR(*edges)[2 * i + 1]);
        if (giant_size != NULL) {VECTOR(*giant_size)[i] = biggest;}
        if (vertex_count != NULL) { VECTOR(*vertex_count)[i] = connected;}
    }

    igraph_vector_int_destroy(&sizes);
    igraph_vector_int_destroy(&links);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bond_percolation
 * \brief The size of the largest component as edges are added to a graph.
 *
 * \experimental
 *
 * Calculates the bond percolation curve, or the size of the largest component as
 * edges are added to the graph in the order given. If both the output and the input
 * are reversed, it is the size of the largest component as edges are removed.
 *
 * \param graph The graph that edges are assumed to be in. Edge directions
 *    are ignored.
 * \param giant_size <code>giant_size[i]</code> will contain the size of the largest
 *    component after having added the edge with index <code>edge_order[i]</code>.
 * \param vertex_count <code>vertex_count[i]</code> will contain the number
 *    of vertices that have at least one incident edge after adding the edge with
 *    index <code>edge_order[i]</code>.
 * \param edge_order The order the edges are added in. Must not contain duplicates.
 *    If \c NULL, a random order will be used.
 * \return Error code.
 *
 * \sa \ref igraph_edgelist_percolation() to specify the edges to be added by
 * their endpoints; \ref igraph_site_percolation() to compute the vertex percolation
 * curve; \ref igraph_connected_components() to find the size of connected components.
 *
 * Time complexity: O(|V| + |E| a(|E|)) where a is the inverse Ackermann function,
 * for all practical purposes it is not above 5.
 */

igraph_error_t igraph_bond_percolation(
        const igraph_t *graph,
        igraph_vector_int_t *output,
        igraph_vector_int_t *vertex_count,
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
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2 * igraph_vector_int_size(p_edge_order));

    IGRAPH_CHECK(igraph_edges(graph, igraph_ess_vector(p_edge_order), &edges));
    IGRAPH_CHECK(igraph_edgelist_percolation(&edges, output, vertex_count));

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
                                     igraph_integer_t *added,
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
        // do not add edges to vertices that have not been added.
        if (VECTOR(*sizes)[VECTOR(*neighbors)[i]] == 0) {
            continue;
        }
        *added += 1;
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
 * Calculates the site percolation curve, or the size of the largest connected
 * component as vertices are added in the order given. If the output is
 * reversed, it is the size of the largest component as vertices are removed
 * in the reverse of the order given. If no vertex order is given, a random
 * one will be used.
 *
 * \param graph The graph that vertices are assumed to be in. Edge directions
 *    are ignored.
 * \param giant_size <code>giant_size[i]</code> will contain the size of the largest
 *    component after having added the vertex with index
 *    <code>vertex_order[i]</code>.
 * \param edge_count <code>edge_count[i]</code> will contain the numer of edges
 *    in the graph having added the vertex with index
 *    <code>vertex_order[i]</code>.
 * \param vertex_order The order the vertices are added in. Must not contain
 *    duplicates. If \c NULL, a random order will be used.
 * \return Error code.
 *
 * \sa \ref igraph_bond_percolation() to compute the edge percolation curve;
 * \ref igraph_connected_components() to find the size of connected components.
 *
 * Time complexity: O(|V| + |E| a(|E|)) where a is the inverse Ackermann function,
 * for all practical purposes it is not above 5.
 */

igraph_error_t igraph_site_percolation(
        const igraph_t *graph,
        igraph_vector_int_t *giant_size,
        igraph_vector_int_t *edge_count,
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
    if (giant_size != NULL) {
        IGRAPH_CHECK(igraph_vector_int_resize(giant_size, vcount));
    }
    if (edge_count != NULL) {
        IGRAPH_CHECK(igraph_vector_int_resize(edge_count, vcount));
    }

    for (igraph_integer_t i = 0; i < vcount; i++) {
        VECTOR(sizes)[i] = 0;
        VECTOR(links)[i] = i;
    }

    igraph_integer_t added = 0;
    igraph_vector_int_t neighbors;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 0);
    for (igraph_integer_t i = 0; i < vcount; i++) {
        IGRAPH_CHECK(percolate_site(graph, &links, &sizes, &biggest, &added, VECTOR(*p_vertex_order)[i], &neighbors));
        if (giant_size != NULL) {
            VECTOR(*giant_size)[i] = biggest;
        }
        if (edge_count != NULL) {
            VECTOR(*edge_count)[i] = added;
        }
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
