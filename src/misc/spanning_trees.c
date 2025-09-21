/*
   igraph library.
   Copyright (C) 2011-2024  The igraph development team <igraph@igraph.org>

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
#include "igraph_bitset.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "core/indheap.h"
#include "core/interruption.h"

static igraph_int_t max_tree_edges(igraph_int_t no_of_nodes, igraph_int_t no_of_edges) {
    if (no_of_nodes == 0) return 0;
    if (no_of_edges < no_of_nodes) return no_of_edges;
    return no_of_nodes - 1;
}


/* Unweighted case -- just a simple BFS */

/**
 * \ingroup structural
 * \function igraph_i_minimum_spanning_tree_unweighted
 * \brief A spanning tree of an unweighted graph.
 *
 * Produces an arbitrary spanning tree of the graph.
 *
 * </para><para>
 * Directed graphs are treated as undirected for this computation.
 *
 * </para><para>
 * If the graph is not connected then a spanning forest is returned.
 * This is a set of spanning trees of each component.
 *
 * \param graph The graph object. Edge directions will be ignored.
 * \param res An initialized vector, the IDs of the edges that constitute
 *        a spanning tree will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanning tree as
 *        a separate graph object.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for temporary data.
 *
 * Time complexity: O(|V|+|E|),
 * |V| is the
 * number of vertices, |E| the number
 * of edges in the graph.
 *
 * \sa \ref igraph_minimum_spanning_tree() for a common interface to all
 * minimum spanning tree algorithms.
 */

static igraph_error_t igraph_i_minimum_spanning_tree_unweighted(
        const igraph_t* graph,
        igraph_vector_int_t* res) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_bitset_t already_added, added_edges;

    igraph_dqueue_int_t q;
    igraph_vector_int_t eids;

    IGRAPH_BITSET_INIT_FINALLY(&added_edges, no_of_edges);
    IGRAPH_BITSET_INIT_FINALLY(&already_added, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    igraph_vector_int_clear(res);
    IGRAPH_CHECK(igraph_vector_int_reserve(res, max_tree_edges(no_of_nodes, no_of_edges)));

    /* Perform a BFS */
    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        if (IGRAPH_BIT_TEST(already_added, i)) {
            continue;
        }

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_BIT_SET(already_added, i);
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        while (! igraph_dqueue_int_empty(&q)) {
            igraph_int_t eids_size;
            igraph_int_t act_node = igraph_dqueue_int_pop(&q);
            IGRAPH_CHECK(igraph_incident(graph, &eids, act_node, IGRAPH_ALL, IGRAPH_LOOPS));
            eids_size = igraph_vector_int_size(&eids);
            for (igraph_int_t j = 0; j < eids_size; j++) {
                igraph_int_t edge = VECTOR(eids)[j];
                if (! IGRAPH_BIT_TEST(added_edges, edge)) {
                    igraph_int_t to = IGRAPH_OTHER(graph, edge, act_node);
                    if (! IGRAPH_BIT_TEST(already_added, to)) {
                        IGRAPH_BIT_SET(already_added, to);
                        IGRAPH_BIT_SET(added_edges, edge);
                        igraph_vector_int_push_back(res, edge); /* reserved */
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, to));
                    }
                }
            }
        }
    }

    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&eids);
    igraph_bitset_destroy(&already_added);
    igraph_bitset_destroy(&added_edges);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/* Prims' algorithm */

/**
 * \ingroup structural
 * \function igraph_i_minimum_spanning_tree_prim
 * \brief A minimum spanning tree of a weighted graph using Prim's method.
 *
 * Finds a spanning tree or spanning forest for which the sum of edge
 * weights is the smallest. This function uses Prim's method for carrying
 * out the computation.
 *
 * </para><para>
 * Directed graphs are treated as undirected for this computation.
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
 * \param res An initialized vector, the IDs of the edges that constitute
 *        a spanning tree will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanning tree as
 *        a separate graph object.
 * \param weights A vector containing the weights of the edges in the order
 *        of edge IDs. Weights must not be NaN.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory.
 *         \c IGRAPH_EINVAL, length of weight vector does not
 *           match number of edges, or NaN in weights.
 *
 * Time complexity: O(|E| log |V|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \sa \ref igraph_minimum_spanning_tree() for a common interface to all
 * minimum spanning tree algorithms.
 *
 * \example examples/simple/igraph_minimum_spanning_tree.c
 */

static igraph_error_t igraph_i_minimum_spanning_tree_prim(
        const igraph_t* graph,
        igraph_vector_int_t* res,
        const igraph_vector_t *weights) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_bitset_t already_added, added_edges;

    igraph_d_indheap_t heap;
    igraph_vector_int_t adj;

    if (weights == NULL) {
        return igraph_i_minimum_spanning_tree_unweighted(graph, res);
    }

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Weight vector length does not match number of edges.", IGRAPH_EINVAL);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
    }

    igraph_vector_int_clear(res);
    IGRAPH_CHECK(igraph_vector_int_reserve(res, max_tree_edges(no_of_nodes, no_of_edges)));

    IGRAPH_BITSET_INIT_FINALLY(&added_edges, no_of_edges);
    IGRAPH_BITSET_INIT_FINALLY(&already_added, no_of_nodes);

    IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
    IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&adj, 0);

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        igraph_int_t adj_size;
        if (IGRAPH_BIT_TEST(already_added, i)) {
            continue;
        }
        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_BIT_SET(already_added, i);
        /* add all edges of the first vertex */
        IGRAPH_CHECK(igraph_incident(graph, &adj, i, IGRAPH_ALL, IGRAPH_LOOPS));
        adj_size = igraph_vector_int_size(&adj);
        for (igraph_int_t j = 0; j < adj_size; j++) {
            igraph_int_t edgeno = VECTOR(adj)[j];
            igraph_int_t neighbor = IGRAPH_OTHER(graph, edgeno, i);
            if (! IGRAPH_BIT_TEST(already_added, neighbor)) {
                IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i, edgeno));
            }
        }

        while (! igraph_d_indheap_empty(&heap)) {
            /* Get minimal edge */
            igraph_int_t from, edge;
            igraph_d_indheap_max_index(&heap, &from, &edge);

            /* Erase it */
            igraph_d_indheap_delete_max(&heap);

            /* Is this edge already included? */
            if (! IGRAPH_BIT_TEST(added_edges, edge)) {
                const igraph_int_t to = IGRAPH_OTHER(graph, edge, from);

                /* Does it point to a visited node? */
                if (! IGRAPH_BIT_TEST(already_added, to)) {
                    IGRAPH_BIT_SET(already_added, to);
                    IGRAPH_BIT_SET(added_edges, edge);
                    igraph_vector_int_push_back(res, edge); /* reserved */
                    /* add all outgoing edges */
                    IGRAPH_CHECK(igraph_incident(graph, &adj, to, IGRAPH_ALL, IGRAPH_LOOPS));
                    adj_size = igraph_vector_int_size(&adj);
                    for (igraph_int_t j = 0; j < adj_size; j++) {
                        const igraph_int_t edgeno = VECTOR(adj)[j];
                        const igraph_int_t neighbor = IGRAPH_OTHER(graph, edgeno, to);
                        if (! IGRAPH_BIT_TEST(already_added, neighbor)) {
                            IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to, edgeno));
                        }
                    }
                } /* for */
            } /* if !already_added */
        } /* while in the same component */
    } /* for all nodes */

    igraph_vector_int_destroy(&adj);
    igraph_d_indheap_destroy(&heap);
    igraph_bitset_destroy(&already_added);
    igraph_bitset_destroy(&added_edges);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/* Kruskal's algorithm */

static igraph_int_t get_comp(igraph_vector_int_t *comp, igraph_int_t i) {
    igraph_int_t k = i;
    for (;;) {
        igraph_int_t next = VECTOR(*comp)[k];
        if (next == k) {
            VECTOR(*comp)[i] = k;
            return k;
        } else {
            k = next;
        }
    }
}

static void merge_comp(igraph_vector_int_t *comp, igraph_int_t i, igraph_int_t j) {
    igraph_int_t ci = get_comp(comp, i);
    igraph_int_t cj = get_comp(comp, j);
    VECTOR(*comp)[ci] = cj;
}

/**
 * \ingroup structural
 * \function igraph_i_minimum_spanning_tree_kruskal
 * \brief A minimum spanning tree of a weighted graph using Kruskal's method.
 *
 * Finds a spanning tree or spanning forest for which the sum of edge
 * weights is the smallest. This function uses Kruskal's method for carrying
 * out the computation.
 *
 * </para><para>
 * Directed graphs are treated as undirected for this computation.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Kruskal, J. B.:
 * On the shortest spanning subtree of a graph and the traveling salesman problem,
 * Proc. Amer. Math. Soc. 7 (1956), 48-50
 * https://doi.org/10.1090%2FS0002-9939-1956-0078686-7
 *
 * \param graph The graph object. Edge directions will be ignored.
 * \param res An initialized vector, the IDs of the edges that constitute
 *        a spanning tree will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanning tree as
 *        a separate graph object.
 * \param weights A vector containing the weights of the edges in the order
 *        of edge IDs. Weights must not be NaN.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory.
 *         \c IGRAPH_EINVAL, length of weight vector does not
 *           match number of edges, or NaN in weights.
 *
 * Time complexity: O(|E| log |E|),
 * |V| is the number of vertices,
 * |E| the number of edges in the
 * graph.
 *
 * \sa \ref igraph_minimum_spanning_tree() for a common interface to all
 * minimum spanning tree algorithms.
 *
 * \example examples/simple/igraph_minimum_spanning_tree.c
 */

static igraph_error_t igraph_i_minimum_spanning_tree_kruskal(
        const igraph_t *graph,
        igraph_vector_int_t *res,
        const igraph_vector_t *weights) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    const igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t idx, comp;
    igraph_int_t tree_edge_count;
    int iter = 0;

    if (weights == NULL) {
        return igraph_i_minimum_spanning_tree_unweighted(graph, res);
    }

    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
        IGRAPH_ERROR("Weight vector length does not match number of edges.", IGRAPH_EINVAL);
    }

    if (igraph_vector_is_any_nan(weights)) {
        IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&idx, no_of_edges);
    IGRAPH_CHECK(igraph_vector_sort_ind(weights, &idx, IGRAPH_ASCENDING));

    IGRAPH_CHECK(igraph_vector_int_init_range(&comp, 0, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &comp);

    igraph_vector_int_clear(res);
    IGRAPH_CHECK(igraph_vector_int_reserve(res, max_tree_edges(no_of_nodes, no_of_edges)));

    tree_edge_count = 0;
    for (igraph_int_t i=0; i < no_of_edges; i++) {
        igraph_int_t edge = VECTOR(idx)[i];
        igraph_int_t u = IGRAPH_FROM(graph, edge);
        igraph_int_t v = IGRAPH_TO(graph, edge);

        igraph_int_t cu = get_comp(&comp, u);
        igraph_int_t cv = get_comp(&comp, v);

        if (cu != cv) {
            merge_comp(&comp, u, v);
            igraph_vector_int_push_back(res, edge); /* reserved */
            tree_edge_count++;
        }

        if (tree_edge_count == no_of_nodes - 1) {
            /* We have enough edges for a tree. */
            break;
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 16);
    }

    igraph_vector_int_destroy(&idx);
    igraph_vector_int_destroy(&comp);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree
 * \brief Calculates a minimum spanning tree of a graph.
 *
 * Finds a minimum weight spanning tree of the graph. If the graph is not
 * connected then its minimum spanning forest is returned, i.e. the set
 * of the minimum spanning trees of each component.
 *
 * </para><para>
 * Directed graphs are treated as undirected for this computation.
 *
 * </para><para>
 * This function is deterministic, i.e. it always returns the same
 * spanning tree. See \ref igraph_random_spanning_tree() for the uniform
 * random sampling of spanning trees of a graph.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Prim, R.C.: Shortest connection networks and some
 * generalizations, Bell System Technical
 * Journal, Vol. 36,
 * 1957, 1389--1401.
 * https://doi.org/10.1002/j.1538-7305.1957.tb01515.x
 *
 * </para><para>
 * Kruskal, J. B.:
 * On the shortest spanning subtree of a graph and the traveling salesman problem,
 * Proc. Amer. Math. Soc. 7 (1956), 48-50
 * https://doi.org/10.1090%2FS0002-9939-1956-0078686-7
 *
 * \param graph The graph object. Edge directions will be ignored.
 * \param res An initialized vector, the IDs of the edges that constitute
 *        a spanning tree will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanning tree as
 *        a separate graph object.
 * \param weights A vector containing the weights of the edges in the order
 *        of edge IDs. Weights must not be NaN. Supply \c NULL to treat all
 *        edges as having the same weight.
 * \param method The type of the algorithm used.
 *        \clist
 *        \cli IGRAPH_MST_AUTOMATIC
 *          tries to select the best performing algorithm for the current graph.
 *        \cli IGRAPH_MST_UNWEIGHTED
 *          ignores edge weights and produces an arbitrary spanning tree.
 *        \cli IGRAPH_MST_PRIM
 *          uses Prim's algorithm.
 *        \cli IGRAPH_MST_KRUSKAL
 *          uses Kruskal's algorithm.
 *        \endclist
 * \return Error code.
 *
 * Time complexity: See the functions implementing the specific algorithms.
 *
 * \sa \ref igraph_random_spanning_tree() to compute a random spanning tree
 *   instead of a minimum one.
 *
 * \example examples/simple/igraph_minimum_spanning_tree.c
 */
igraph_error_t igraph_minimum_spanning_tree(
    const igraph_t *graph, igraph_vector_int_t *res,
    const igraph_vector_t *weights, igraph_mst_algorithm_t method)
{
    if (method == IGRAPH_MST_AUTOMATIC) {
        /* For now we use igraph_minimum_spanning_tree_kruskal() unconditionally
         * for the weighted case; see benchmarks. */
        if (weights == NULL) {
            method = IGRAPH_MST_UNWEIGHTED;
        } else {
            method = IGRAPH_MST_KRUSKAL;
        }
    }
    switch (method) {
    case IGRAPH_MST_UNWEIGHTED:
        return igraph_i_minimum_spanning_tree_unweighted(graph, res);
    case IGRAPH_MST_PRIM:
        return igraph_i_minimum_spanning_tree_prim(graph, res, weights);
    case IGRAPH_MST_KRUSKAL:
        return igraph_i_minimum_spanning_tree_kruskal(graph, res, weights);
    default:
        IGRAPH_ERROR("Invalid method for minimum spanning tree.", IGRAPH_EINVAL);
    }
}


/* igraph_random_spanning_tree */

/* Loop-erased random walk (LERW) implementation.
 * res must be an initialized vector. The edge IDs of the spanning tree
 * will be added to the end of it. res will not be cleared before doing this.
 *
 * The walk is started from vertex start. comp_size must be the size of the connected
 * component containing start.
 */
static igraph_error_t igraph_i_lerw(const igraph_t *graph, igraph_vector_int_t *res, igraph_int_t start,
                         igraph_int_t comp_size, igraph_vector_bool_t *visited, const igraph_inclist_t *il) {
    igraph_int_t visited_count;

    IGRAPH_CHECK(igraph_vector_int_reserve(res, igraph_vector_int_size(res) + comp_size - 1));

    VECTOR(*visited)[start] = true;
    visited_count = 1;

    while (visited_count < comp_size) {
        igraph_int_t degree, edge;
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
igraph_error_t igraph_random_spanning_tree(const igraph_t *graph, igraph_vector_int_t *res, igraph_int_t vid) {
    igraph_inclist_t il;
    igraph_vector_bool_t visited;
    igraph_int_t vcount = igraph_vcount(graph);

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
        igraph_int_t comp_count;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&membership, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&csize, 0);

        IGRAPH_CHECK(igraph_connected_components(graph, &membership, &csize, &comp_count, IGRAPH_WEAK));

        /* for each component ... */
        for (igraph_int_t i = 0; i < comp_count; ++i) {
            /* ... find a vertex to start the LERW from */
            igraph_int_t j = 0;
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
        igraph_int_t comp_size;

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
