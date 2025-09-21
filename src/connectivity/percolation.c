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

#include "igraph_components.h"

#include "igraph_bitset.h"
#include "igraph_constants.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "core/interruption.h"

/**
 * \function percolate_edge
 * \brief Percolates a single edge.
 *
 * \param links Vector representing parents.
 * \param sizes sizes[i] is the number of children of links[i]
 * \param biggest The biggest value in sizes, is updated if a bigger cluster is created.
 * \param a A vertex incident to the edge.
 * \param b The other vertex incident to the edge.
 */

static void percolate_edge(igraph_vector_int_t *links,
                           igraph_vector_int_t *sizes,
                           igraph_int_t *biggest,
                           igraph_int_t a,
                           igraph_int_t b) {

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
    igraph_int_t parent, child;
    if (VECTOR(*sizes)[a] < VECTOR(*sizes)[b]) {
        parent = b;
        child = a;
    } else {
        parent = a;
        child = b;
    }

    VECTOR(*links)[child] = parent;
    VECTOR(*sizes)[parent] += VECTOR(*sizes)[child];

    // If made new biggest component, update biggest
    if (VECTOR(*sizes)[parent] >= *biggest) {
        *biggest = VECTOR(*sizes)[parent];
    }
}

/**
 * \function igraph_edgelist_percolation
 * \brief The size of the largest component as vertex pairs are connected.
 *
 * \experimental
 *
 * Calculates the size of the largest connected component as edges are added
 * to a graph in the given order. This function differs from
 * \ref igraph_bond_percolation() in that it take a list of vertex pairs as input.
 *
 * \param edges Vector of edges, where the i-th edge has endpoints
 *    <code>edges[2i]</code> and <code>edges[2i+1]</code>.
 * \param giant_size <code>giant_size[i]</code> will contain the size of the
 *    largest connected component after edge \c i is added.
 * \param vertex_count <code>vertex_count[i]</code> will contain the number of
 *    vertices with at least one edge after edge \c i is added.
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

    igraph_int_t biggest = 1;
    igraph_int_t vertices_added = 0;
    igraph_int_t lower, upper;
    int iter = 0;

    igraph_int_t ecount = igraph_vector_int_size(edges);

    if (ecount % 2 == 1) {
        IGRAPH_ERROR("Invalid edge list, odd number of elements.", IGRAPH_EINVAL);
    }
    ecount = ecount / 2;

    if (giant_size != NULL) {
        IGRAPH_CHECK(igraph_vector_int_resize(giant_size, ecount));
    }
    if (vertex_count != NULL) {
        IGRAPH_CHECK(igraph_vector_int_resize(vertex_count, ecount));
    }

    // Handle edge case of no edges.
    if (ecount == 0) {
         return IGRAPH_SUCCESS;
    }

    igraph_vector_int_minmax(edges, &lower, &upper);

    if (lower < 0) {
        IGRAPH_ERROR("Invalid vertex ID.", IGRAPH_EINVVID);
    }

    const igraph_int_t vcount = upper + 1;

    igraph_vector_int_t sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, vcount);

    igraph_vector_int_t links;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&links, vcount);

    for (igraph_int_t i = 0; i < vcount; i++) {
        VECTOR(sizes)[i] = -1;
        VECTOR(links)[i] =  i;
    }

    for (igraph_int_t i = 0; i < ecount; i++) {
        const igraph_int_t from = VECTOR(*edges)[2*i];
        const igraph_int_t to   = VECTOR(*edges)[2*i + 1];
        if (VECTOR(sizes)[from] == -1) {
            vertices_added++;
            VECTOR(sizes)[from] = 1;
        }
        if (VECTOR(sizes)[to] == -1) {
            vertices_added++;
            VECTOR(sizes)[to] = 1;
        }
        percolate_edge(&links, &sizes, &biggest, from, to);
        if (giant_size != NULL) {
            VECTOR(*giant_size)[i] = biggest;
        }
        if (vertex_count != NULL) {
            VECTOR(*vertex_count)[i] = vertices_added;
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }

    igraph_vector_int_destroy(&links);
    igraph_vector_int_destroy(&sizes);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bond_percolation
 * \brief The size of the largest component as edges are added to a graph.
 *
 * \experimental
 *
 * Calculates the bond percolation curve, i.e. the size of the largest connected
 * component as edges are added to the graph in the order given. If both
 * \p giant_size and \p edge_order are reversed, it is the size of the largest
 * component as edges are removed from the graph. If no edge order is given,
 * a random one will be used.
 *
 * \param graph The graph that edges are assumed to be in. Edge directions
 *    are ignored.
 * \param giant_size <code>giant_size[i]</code> will contain the size of the
 *    largest component after having added the edge with index
 *    <code>edge_order[i]</code>.
 * \param vertex_count <code>vertex_count[i]</code> will contain the number
 *    of vertices that have at least one incident edge after adding the edge
 *    with index <code>edge_order[i]</code>.
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
        igraph_vector_int_t *giant_size,
        igraph_vector_int_t *vertex_count,
        const igraph_vector_int_t *edge_order) {

    const igraph_vector_int_t *p_edge_order;
    igraph_vector_int_t i_edge_order;
    igraph_vector_int_t edges;

    // Use a random edge order when no edge order was given
    if (edge_order == NULL) {
        IGRAPH_CHECK(igraph_vector_int_init_range(&i_edge_order, 0, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &i_edge_order);
        igraph_vector_int_shuffle(&i_edge_order);
        p_edge_order = &i_edge_order;
    } else {
        // Verify that there are no duplicates.
        const igraph_int_t no_of_added_edges = igraph_vector_int_size(edge_order);
        igraph_bitset_t present_edges;

        IGRAPH_BITSET_INIT_FINALLY(&present_edges, no_of_added_edges);

        for (igraph_int_t i = 0; i < no_of_added_edges; i++) {
            if (IGRAPH_BIT_TEST(present_edges, VECTOR(*edge_order)[i])) {
                IGRAPH_ERROR("Duplicate edges in edge order vector.", IGRAPH_EINVAL);
            }
            IGRAPH_BIT_SET(present_edges, VECTOR(*edge_order)[i]);
        }
        igraph_bitset_destroy(&present_edges);
        IGRAPH_FINALLY_CLEAN(1);
        p_edge_order = edge_order;
    }

    // Initialize edge list. igraph_edges() will validate edge IDs.
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2 * igraph_vector_int_size(p_edge_order));
    IGRAPH_CHECK(igraph_edges(graph, igraph_ess_vector(p_edge_order), &edges, /* bycol = */ 0));

    // Defer to igraph_edgelist_percolation()
    IGRAPH_CHECK(igraph_edgelist_percolation(&edges, giant_size, vertex_count));

    // Cleanup
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
                                     igraph_int_t *biggest,
                                     igraph_int_t *edges_added,
                                     igraph_int_t vertex,
                                     igraph_vector_int_t *neighbors) {

    if (VECTOR(*sizes)[vertex] != 0) {
        IGRAPH_ERROR("Duplicate vertices in vertex order vector.", IGRAPH_EINVAL);
    }

    VECTOR(*sizes)[vertex] = 1;

    IGRAPH_CHECK(igraph_neighbors(graph, neighbors, vertex, IGRAPH_ALL, IGRAPH_LOOPS, IGRAPH_MULTIPLE));

    igraph_int_t neighbor_count = igraph_vector_int_size(neighbors);
    for (igraph_int_t i = 0; i < neighbor_count; i++) {
        // Do not add edges to vertices that have not been added.
        if (VECTOR(*sizes)[VECTOR(*neighbors)[i]] == 0) {
            continue;
        }
        *edges_added += 1;
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
 * Calculates the site percolation curve, i.e. the size of the largest connected
 * component as vertices are added in the given order. If both \p giant_size
 * and \p vertex_order are reversed, it is the size of the largest component
 * as vertices are removed from the graph. If no vertex order is given, a random
 * one will be used.
 *
 * \param graph The graph that vertices are assumed to be in. Edge directions
 *    are ignored.
 * \param giant_size <code>giant_size[i]</code> will contain the size of the
 *    largest component after having added the vertex with index
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

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_vector_int_t *p_vertex_order;
    igraph_vector_int_t i_vertex_order;
    int iter = 0;

    // Use a random vertex order when no vertex order was given
    if (vertex_order == NULL) {
        IGRAPH_CHECK(igraph_vector_int_init_range(&i_vertex_order, 0, vcount));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &i_vertex_order);
        igraph_vector_int_shuffle(&i_vertex_order);
        p_vertex_order = &i_vertex_order;
    } else {
        p_vertex_order = vertex_order;
    }

    // Initialize variables
    igraph_int_t number_percolated = igraph_vector_int_size(p_vertex_order);
    igraph_int_t biggest = 1; // largest component size so far
    igraph_int_t edges_added = 0; // no. of edges added so far

    igraph_vector_int_t sizes;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, vcount);

    igraph_vector_int_t links;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&links, vcount);

    for (igraph_int_t i = 0; i < vcount; i++) {
        VECTOR(sizes)[i] = 0;
        VECTOR(links)[i] = i;
    }

    igraph_vector_int_t neighbors;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 0);

    if (giant_size != NULL) {
        IGRAPH_CHECK(igraph_vector_int_resize(giant_size, number_percolated));
    }
    if (edge_count != NULL) {
        IGRAPH_CHECK(igraph_vector_int_resize(edge_count, number_percolated));
    }

    // Percolation
    for (igraph_int_t i = 0; i < number_percolated; i++) {
        const igraph_int_t vid = VECTOR(*p_vertex_order)[i];
        if (vid < 0 || vid >= vcount) {
            IGRAPH_ERROR("Invalid vertex ID.", IGRAPH_EINVVID);
        }
        IGRAPH_CHECK(percolate_site(graph, &links, &sizes,
                                    &biggest, &edges_added, vid, &neighbors));
        if (giant_size != NULL) {
            VECTOR(*giant_size)[i] = biggest;
        }
        if (edge_count != NULL) {
            VECTOR(*edge_count)[i] = edges_added;
        }
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }

    // Cleanup
    igraph_vector_int_destroy(&neighbors);
    igraph_vector_int_destroy(&links);
    igraph_vector_int_destroy(&sizes);
    IGRAPH_FINALLY_CLEAN(3);

    if (vertex_order == NULL) {
        igraph_vector_int_destroy(&i_vertex_order);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
