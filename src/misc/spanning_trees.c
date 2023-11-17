/*
   IGraph library.
   Copyright (C) 2011-2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "core/indheap.h"
#include "core/interruption.h"

static igraph_error_t igraph_i_minimum_spanning_tree_unweighted(
    const igraph_t *graph, igraph_vector_int_t *result);
static igraph_error_t igraph_i_minimum_spanning_tree_prim(
    const igraph_t *graph, igraph_vector_int_t *result, const igraph_vector_t *weights);

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree
 * \brief Calculates one minimum spanning tree of a graph.
 *
 * Finds a spanning tree of the graph. If the graph is not connected
 * then its minimum spanning forest is returned. This is the set of the
 * minimum spanning trees of each component.
 *
 * </para><para>
 * Directed graphs are considered as undirected for this computation.
 *
 * </para><para>
 * This function is deterministic, i.e. it always returns the same
 * spanning tree. See \ref igraph_random_spanning_tree() for the uniform
 * random sampling of spanning trees of a graph.
 *
 * \param graph The graph object.
 * \param res An initialized vector, the IDs of the edges that constitute
 *        a spanning tree will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanning tree as
 *        a separate graph object.
 * \param weights A vector containing the weights of the edges
 *        in the same order as the simple edge iterator visits them
 *        (i.e. in increasing order of edge IDs).
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V|+|E|) for the unweighted case, O(|E| log |V|)
 * for the weighted case. |V| is the number of vertices, |E| the
 * number of edges in the graph.
 *
 * \sa \ref igraph_minimum_spanning_tree_unweighted() and
 *     \ref igraph_minimum_spanning_tree_prim() if you only need the
 *     tree as a separate graph object.
 *
 * \example examples/simple/igraph_minimum_spanning_tree.c
 */
igraph_error_t igraph_minimum_spanning_tree(
    const igraph_t *graph, igraph_vector_int_t *res, const igraph_vector_t *weights
) {
    if (weights == NULL) {
        IGRAPH_CHECK(igraph_i_minimum_spanning_tree_unweighted(graph, res));
    } else {
        IGRAPH_CHECK(igraph_i_minimum_spanning_tree_prim(graph, res, weights));
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_unweighted
 * \brief Calculates one minimum spanning tree of an unweighted graph.
 *
 * If the graph has more minimum spanning trees (this is always the
 * case, except if it is a forest) this implementation returns only
 * the same one.
 *
 * </para><para>
 * Directed graphs are considered as undirected for this computation.
 *
 * </para><para>
 * If the graph is not connected then its minimum spanning forest is
 * returned. This is the set of the minimum spanning trees of each
 * component.
 *
 * \param graph The graph object. Edge directions will be ignored.
 * \param mst The minimum spanning tree, another graph object. Do
 *        \em not initialize this object before passing it to
 *        this function, but be sure to call \ref igraph_destroy() on it if
 *        you don't need it any more.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the
 * number of vertices, |E| the number
 * of edges in the graph.
 *
 * \sa \ref igraph_minimum_spanning_tree_prim() for weighted graphs,
 *     \ref igraph_minimum_spanning_tree() if you need the IDs of the
 *     edges that constitute the spanning tree.
 */

igraph_error_t igraph_minimum_spanning_tree_unweighted(const igraph_t *graph,
        igraph_t *mst) {
    igraph_vector_int_t edges;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_nodes > 0 ? no_of_nodes - 1 : 0);
    IGRAPH_CHECK(igraph_i_minimum_spanning_tree_unweighted(graph, &edges));
    IGRAPH_CHECK(igraph_subgraph_from_edges(
        graph, mst, igraph_ess_vector(&edges), /* delete_vertices = */ false));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_prim
 * \brief Calculates one minimum spanning tree of a weighted graph.
 *
 * Finds a spanning tree or spanning forest for which the sum of edge
 * weights is the smallest. This function uses Prim's method for carrying
 * out the computation.
 *
 * </para><para>
 * Directed graphs are considered as undirected for this computation.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Prim, R.C.: Shortest connection networks and some
 * generalizations, Bell System Technical
 * Journal, Vol. 36,
 * 1957, 1389--1401.
 * https://doi.org/10.1002/j.1538-7305.1957.tb01515.x
 *
 * \param graph The graph object. Edge directions will be ignored.
 * \param mst The result of the computation, a graph object containing
 *        the minimum spanning tree of the graph.
 *        Do \em not initialize this object before passing it to
 *        this function, but be sure to call \ref igraph_destroy() on it if
 *        you don't need it any more.
 * \param weights A vector containing the weights of the edges
 *        in the same order as the simple edge iterator visits them
 *        (i.e. in increasing order of edge IDs).
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory.
 *         \c IGRAPH_EINVAL, length of weight vector does not
 *           match number of edges.
 *
 * Time complexity: O(|E| log |V|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \sa \ref igraph_minimum_spanning_tree_unweighted() for unweighted graphs,
 *     \ref igraph_minimum_spanning_tree() if you need the IDs of the
 *     edges that constitute the spanning tree.
 *
 * \example examples/simple/igraph_minimum_spanning_tree.c
 */

igraph_error_t igraph_minimum_spanning_tree_prim(const igraph_t *graph, igraph_t *mst,
                                      const igraph_vector_t *weights) {
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, igraph_vcount(graph) - 1);
    IGRAPH_CHECK(igraph_i_minimum_spanning_tree_prim(graph, &edges, weights));
    IGRAPH_CHECK(igraph_subgraph_from_edges(
        graph, mst, igraph_ess_vector(&edges), /* delete_vertices = */ false));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_i_minimum_spanning_tree_unweighted(const igraph_t* graph, igraph_vector_int_t* res) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    bool *already_added, *added_edges;

    igraph_dqueue_int_t q;
    igraph_vector_int_t eids;

    igraph_vector_int_clear(res);

    added_edges = IGRAPH_CALLOC(no_of_edges, bool);
    IGRAPH_CHECK_OOM(added_edges, "Insufficient memory for unweighted spanning tree.");
    IGRAPH_FINALLY(igraph_free, added_edges);

    already_added = IGRAPH_CALLOC(no_of_nodes, bool);
    IGRAPH_CHECK_OOM(already_added, "Insufficient memory for unweighted spanning tree.");
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    /* Perform a BFS */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        if (already_added[i]) {
            continue;
        }

        IGRAPH_ALLOW_INTERRUPTION();

        already_added[i] = true;
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        while (! igraph_dqueue_int_empty(&q)) {
            igraph_integer_t eids_size;
            igraph_integer_t act_node = igraph_dqueue_int_pop(&q);
            IGRAPH_CHECK(igraph_incident(graph, &eids, act_node,
                                         IGRAPH_ALL));
            eids_size = igraph_vector_int_size(&eids);
            for (igraph_integer_t j = 0; j < eids_size; j++) {
                igraph_integer_t edge = VECTOR(eids)[j];
                if (! added_edges[edge]) {
                    igraph_integer_t to = IGRAPH_OTHER(graph, edge, act_node);
                    if (! already_added[to]) {
                        already_added[to] = true;
                        added_edges[edge] = true;
                        IGRAPH_CHECK(igraph_vector_int_push_back(res, edge));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, to));
                    }
                }
            }
        }
    }

    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&eids);
    IGRAPH_FREE(already_added);
    IGRAPH_FREE(added_edges);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_minimum_spanning_tree_prim(
        const igraph_t* graph, igraph_vector_int_t* res, const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    bool *already_added, *added_edges;

    igraph_d_indheap_t heap;
    const igraph_neimode_t mode = IGRAPH_ALL;

    igraph_vector_int_t adj;

    igraph_vector_int_clear(res);

    if (weights == NULL) {
        return igraph_i_minimum_spanning_tree_unweighted(graph, res);
    }

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Weight vector length does not match number of edges.", IGRAPH_EINVAL);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weigths must not contain NaN values.", IGRAPH_EINVAL);
    }

    added_edges = IGRAPH_CALLOC(no_of_edges, bool);
    IGRAPH_CHECK_OOM(added_edges, "Insufficient memory for minimum spanning tree calculation.");
    IGRAPH_FINALLY(igraph_free, added_edges);

    already_added = IGRAPH_CALLOC(no_of_nodes, bool);
    IGRAPH_CHECK_OOM(already_added, "Insufficient memory for minimum spanning tree calculation.");
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
    IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&adj, 0);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t adj_size;
        if (already_added[i]) {
            continue;
        }
        IGRAPH_ALLOW_INTERRUPTION();

        already_added[i] = true;
        /* add all edges of the first vertex */
        IGRAPH_CHECK(igraph_incident(graph, &adj, i, mode));
        adj_size = igraph_vector_int_size(&adj);
        for (igraph_integer_t j = 0; j < adj_size; j++) {
            igraph_integer_t edgeno = VECTOR(adj)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(graph, edgeno, i);
            if (! already_added[neighbor]) {
                IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i, edgeno));
            }
        }

        while (! igraph_d_indheap_empty(&heap)) {
            /* Get minimal edge */
            igraph_integer_t from, edge;
            igraph_d_indheap_max_index(&heap, &from, &edge);

            /* Erase it */
            igraph_d_indheap_delete_max(&heap);

            /* Is this edge already included? */
            if (! added_edges[edge]) {
                igraph_integer_t to = IGRAPH_OTHER(graph, edge, from);

                /* Does it point to a visited node? */
                if (! already_added[to]) {
                    already_added[to] = true;
                    added_edges[edge] = true;
                    IGRAPH_CHECK(igraph_vector_int_push_back(res, edge));
                    /* add all outgoing edges */
                    IGRAPH_CHECK(igraph_incident(graph, &adj, to, mode));
                    adj_size = igraph_vector_int_size(&adj);
                    for (igraph_integer_t j = 0; j < adj_size; j++) {
                        igraph_integer_t edgeno = VECTOR(adj)[j];
                        igraph_integer_t neighbor = IGRAPH_OTHER(graph, edgeno, to);
                        if (! already_added[neighbor]) {
                            IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to, edgeno));
                        }
                    }
                } /* for */
            } /* if !already_added */
        } /* while in the same component */
    } /* for all nodes */

    igraph_d_indheap_destroy(&heap);
    IGRAPH_FREE(already_added);
    igraph_vector_int_destroy(&adj);
    IGRAPH_FREE(added_edges);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/* igraph_random_spanning_tree */

/* Loop-erased random walk (LERW) implementation.
 * res must be an initialized vector. The edge IDs of the spanning tree
 * will be added to the end of it. res will not be cleared before doing this.
 *
 * The walk is started from vertex start. comp_size must be the size of the connected
 * component containing start.
 */
static igraph_error_t igraph_i_lerw(const igraph_t *graph, igraph_vector_int_t *res, igraph_integer_t start,
                         igraph_integer_t comp_size, igraph_vector_bool_t *visited, const igraph_inclist_t *il) {
    igraph_integer_t visited_count;

    IGRAPH_CHECK(igraph_vector_int_reserve(res, igraph_vector_int_size(res) + comp_size - 1));

    VECTOR(*visited)[start] = true;
    visited_count = 1;

    RNG_BEGIN();

    while (visited_count < comp_size) {
        igraph_integer_t degree, edge;
        igraph_vector_int_t *edges;

        edges = igraph_inclist_get(il, start);

        /* choose a random edge */
        degree = igraph_vector_int_size(edges);
        edge = VECTOR(*edges)[ RNG_INTEGER(0, degree - 1) ];

        /* set 'start' to the next vertex */
        start = IGRAPH_OTHER(graph, edge, start);

        /* if the next vertex hasn't been visited yet, register the edge we just traversed */
        if (! VECTOR(*visited)[start]) {
            IGRAPH_CHECK(igraph_vector_int_push_back(res, edge));
            VECTOR(*visited)[start] = true;
            visited_count++;
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_random_spanning_tree
 * \brief Uniformly samples the spanning trees of a graph.
 *
 * Performs a loop-erased random walk on the graph to uniformly sample
 * its spanning trees. Edge directions are ignored.
 * </para><para>
 *
 * Multi-graphs are supported, and edge multiplicities will affect the sampling
 * frequency. For example, consider the 3-cycle graph <code>1=2-3-1</code>, with two edges
 * between vertices 1 and 2. Due to these parallel edges, the trees <code>1-2-3</code>
 * and <code>3-1-2</code> will be sampled with multiplicity 2, while the tree
 * <code>2-3-1</code> will be sampled with multiplicity 1.
 *
 * \param graph The input graph. Edge directions are ignored.
 * \param res An initialized vector, the IDs of the edges that constitute
 *        a spanning tree will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanning tree as
 *        a separate graph object.
 * \param vid This parameter is relevant if the graph is not connected.
 *        If negative, a random spanning forest of all components will be
 *        generated. Otherwise, it should be the ID of a vertex. A random
 *        spanning tree of the component containing the vertex will be
 *        generated.
 *
 * \return Error code.
 *
 * \sa \ref igraph_minimum_spanning_tree(), \ref igraph_random_walk()
 *
 */
igraph_error_t igraph_random_spanning_tree(const igraph_t *graph, igraph_vector_int_t *res, igraph_integer_t vid) {
    igraph_inclist_t il;
    igraph_vector_bool_t visited;
    igraph_integer_t vcount = igraph_vcount(graph);

    if (vid >= vcount) {
        IGRAPH_ERROR("Invalid vertex ID given for random spanning tree.", IGRAPH_EINVVID);
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &il, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &il);

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, vcount));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    igraph_vector_int_clear(res);

    if (vid < 0) { /* generate random spanning forest: consider each component separately */
        igraph_vector_int_t membership, csize;
        igraph_integer_t comp_count;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&csize, 0);

        IGRAPH_CHECK(igraph_connected_components(graph, &membership, &csize, &comp_count, IGRAPH_WEAK));

        /* for each component ... */
        for (igraph_integer_t i = 0; i < comp_count; ++i) {
            /* ... find a vertex to start the LERW from */
            igraph_integer_t j = 0;
            while (VECTOR(membership)[j] != i) {
                ++j;
            }

            IGRAPH_CHECK(igraph_i_lerw(graph, res, j, VECTOR(csize)[i], &visited, &il));
        }

        igraph_vector_int_destroy(&membership);
        igraph_vector_int_destroy(&csize);
        IGRAPH_FINALLY_CLEAN(2);
    } else { /* consider the component containing vid */
        igraph_vector_int_t comp_vertices;
        igraph_integer_t comp_size;

        /* we measure the size of the component */
        IGRAPH_VECTOR_INT_INIT_FINALLY(&comp_vertices, 0);
        IGRAPH_CHECK(igraph_subcomponent(graph, &comp_vertices, vid, IGRAPH_ALL));
        comp_size = igraph_vector_int_size(&comp_vertices);
        igraph_vector_int_destroy(&comp_vertices);
        IGRAPH_FINALLY_CLEAN(1);

        IGRAPH_CHECK(igraph_i_lerw(graph, res, vid, comp_size, &visited, &il));
    }

    igraph_vector_bool_destroy(&visited);
    igraph_inclist_destroy(&il);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
