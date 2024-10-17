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
#include "igraph_bitset.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_operators.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "core/indheap.h"
#include "core/interruption.h"

static igraph_error_t igraph_i_minimum_spanning_tree_unweighted(
    const igraph_t *graph, igraph_vector_int_t *result);
static igraph_error_t igraph_i_minimum_spanning_tree_prim(
    const igraph_t *graph, igraph_vector_int_t *result, const igraph_vector_t *weights);

typedef enum {
  IGRAPH_MST_UNWEIGHTED,
  IGRAPH_MST_PRIM,
  IGRAPH_MST_KRUSKAL
} igraph_mst_algorithm_t;

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
 * \sa \ref igraph_minimum_spanning_tree_unweighted() for unweighted graphs,
 *     \ref igraph_minimum_spanning_tree_prim() for weighted graphs, and
 *     \ref igraph_random_spanning_tree_kruskal() for weighted graphs.
 *
 * \example examples/simple/igraph_minimum_spanning_tree.c
 */

igraph_error_t igraph_minimum_spanning_tree(const igraph_t *graph,
                                            igraph_vector_int_t *res,
                                            const igraph_vector_t *weights,
                                            igraph_mst_algorithm_t algorithm) {
  if (weights == NULL) {
    // Use the unweighted MST algorithm if no weights are provided
    IGRAPH_CHECK(igraph_i_minimum_spanning_tree_unweighted(graph, res));
  } else {
    // Select the appropriate algorithm based on the parameter
    switch (algorithm) {
      case IGRAPH_MST_PRIM:
        IGRAPH_CHECK(igraph_i_minimum_spanning_tree_prim(graph, res, weights));
        break;
      case IGRAPH_MST_KRUSKAL:
        IGRAPH_CHECK(
            igraph_i_minimum_spanning_tree_kruskal(graph, res, weights));
        break;
      default:
        return IGRAPH_EINVAL;  // Invalid algorithm specified
    }
  }
  return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_unweighted
 * \brief Calculates one minimum spanning tree of an unweighted graph.
 *
 * \deprecated-by igraph_minimum_spanning_tree 0.10.14
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
 * \deprecated-by igraph_minimum_spanning_tree 0.10.14
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
    igraph_bitset_t already_added, added_edges;

    igraph_dqueue_int_t q;
    igraph_vector_int_t eids;

    igraph_vector_int_clear(res);

    IGRAPH_BITSET_INIT_FINALLY(&added_edges, no_of_edges);
    IGRAPH_BITSET_INIT_FINALLY(&already_added, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    /* Perform a BFS */
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        if (IGRAPH_BIT_TEST(already_added, i)) {
            continue;
        }

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_BIT_SET(already_added, i);
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        while (! igraph_dqueue_int_empty(&q)) {
            igraph_integer_t eids_size;
            igraph_integer_t act_node = igraph_dqueue_int_pop(&q);
            IGRAPH_CHECK(igraph_incident(graph, &eids, act_node,
                                         IGRAPH_ALL));
            eids_size = igraph_vector_int_size(&eids);
            for (igraph_integer_t j = 0; j < eids_size; j++) {
                igraph_integer_t edge = VECTOR(eids)[j];
                if (! IGRAPH_BIT_TEST(added_edges, edge)) {
                    igraph_integer_t to = IGRAPH_OTHER(graph, edge, act_node);
                    if (! IGRAPH_BIT_TEST(already_added, to)) {
                        IGRAPH_BIT_SET(already_added, to);
                        IGRAPH_BIT_SET(added_edges, edge);
                        IGRAPH_CHECK(igraph_vector_int_push_back(res, edge));
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

static igraph_error_t igraph_i_minimum_spanning_tree_prim(
        const igraph_t* graph, igraph_vector_int_t* res, const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bitset_t already_added, added_edges;

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

    IGRAPH_BITSET_INIT_FINALLY(&added_edges, no_of_edges);
    IGRAPH_BITSET_INIT_FINALLY(&already_added, no_of_nodes);

    IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
    IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&adj, 0);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t adj_size;
        if (IGRAPH_BIT_TEST(already_added, i)) {
            continue;
        }
        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_BIT_SET(already_added, i);
        /* add all edges of the first vertex */
        IGRAPH_CHECK(igraph_incident(graph, &adj, i, mode));
        adj_size = igraph_vector_int_size(&adj);
        for (igraph_integer_t j = 0; j < adj_size; j++) {
            igraph_integer_t edgeno = VECTOR(adj)[j];
            igraph_integer_t neighbor = IGRAPH_OTHER(graph, edgeno, i);
            if (! IGRAPH_BIT_TEST(already_added, neighbor)) {
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
            if (! IGRAPH_BIT_TEST(added_edges, edge)) {
                igraph_integer_t to = IGRAPH_OTHER(graph, edge, from);

                /* Does it point to a visited node? */
                if (! IGRAPH_BIT_TEST(already_added, to)) {
                    IGRAPH_BIT_SET(already_added, to);
                    IGRAPH_BIT_SET(added_edges, edge);
                    IGRAPH_CHECK(igraph_vector_int_push_back(res, edge));
                    /* add all outgoing edges */
                    IGRAPH_CHECK(igraph_incident(graph, &adj, to, mode));
                    adj_size = igraph_vector_int_size(&adj);
                    for (igraph_integer_t j = 0; j < adj_size; j++) {
                        igraph_integer_t edgeno = VECTOR(adj)[j];
                        igraph_integer_t neighbor = IGRAPH_OTHER(graph, edgeno, to);
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

typedef struct {
  igraph_vector_int_t parent;
  igraph_vector_int_t rank;
} igraph_union_find_t;

// Initialize Union-Find structure
static igraph_error_t igraph_union_find_init(igraph_union_find_t *uf,
                                             igraph_integer_t n) {
  IGRAPH_CHECK(igraph_vector_int_init(&uf->parent, n));
  IGRAPH_CHECK(igraph_vector_int_init(&uf->rank, n));
  for (igraph_integer_t i = 0; i < n; i++) {
    VECTOR(uf->parent)[i] = i;  // each node is its own parent initially
    VECTOR(uf->rank)[i] = 0;    // rank is initially zero
  }
  return IGRAPH_SUCCESS;
}

// Find the root of the set containing element i with path compression
static igraph_integer_t igraph_union_find_find(igraph_union_find_t *uf,
                                               igraph_integer_t i) {
  if (VECTOR(uf->parent)[i] != i) {
    VECTOR(uf->parent)
    [i] =
        igraph_union_find_find(uf, VECTOR(uf->parent)[i]);  // path compression
  }
  return VECTOR(uf->parent)[i];
}

// Union of two sets containing elements x and y, by rank
static void igraph_union_find_union(igraph_union_find_t *uf, igraph_integer_t x,
                                    igraph_integer_t y) {
  igraph_integer_t root_x = igraph_union_find_find(uf, x);
  igraph_integer_t root_y = igraph_union_find_find(uf, y);

  if (root_x != root_y) {
    // Union by rank
    if (VECTOR(uf->rank)[root_x] > VECTOR(uf->rank)[root_y]) {
      VECTOR(uf->parent)[root_y] = root_x;
    } else if (VECTOR(uf->rank)[root_x] < VECTOR(uf->rank)[root_y]) {
      VECTOR(uf->parent)[root_x] = root_y;
    } else {
      VECTOR(uf->parent)[root_y] = root_x;
      VECTOR(uf->rank)[root_x]++;
    }
  }
}

// Destroy Union-Find structure
static void igraph_union_find_destroy(igraph_union_find_t *uf) {
  igraph_vector_int_destroy(&uf->parent);
  igraph_vector_int_destroy(&uf->rank);
}

/**
 * \ingroup structural
 * \function igraph_minimum_spanning_tree_kruskal
 * \brief Calculates one minimum spanning tree of a weighted graph.
 *
 * \deprecated-by igraph_minimum_spanning_tree 0.10.14
 *
 * Finds a spanning tree or spanning forest for which the sum of edge
 * weights is the smallest. This function uses Kruskal's method for carrying
 * out the computation.
 *
 * </para><para>
 * Directed graphs are considered as undirected for this computation.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Kruskal, J.B.: On the shortest spanning subtree of a graph and the traveling
 * salesman problem, American Mathematical Society
 * https://doi.org/10.1090/S0002-9939-1956-0078686-7
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

igraph_error_t igraph_minimum_spanning_tree_kruskal(
    const igraph_t *graph, igraph_t *mst, const igraph_vector_t *weights) {
  igraph_vector_int_t edges;

  IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, igraph_vcount(graph) - 1);
  IGRAPH_CHECK(igraph_i_minimum_spanning_tree_kruskal(graph, &edges, weights));
  IGRAPH_CHECK(igraph_subgraph_from_edges(graph, mst, igraph_ess_vector(&edges),
                                          /* delete_vertices = */ false));
  igraph_vector_int_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_minimum_spanning_tree_kruskal(
    const igraph_t *graph, igraph_vector_int_t *res,
    const igraph_vector_t *weights) {
  igraph_integer_t no_of_edges = igraph_ecount(graph);
  igraph_integer_t no_of_nodes = igraph_vcount(graph);

  // Edge list and corresponding weights
  igraph_vector_t sorted_edges;
  IGRAPH_VECTOR_INIT_FINALLY(&sorted_edges, no_of_edges);

  igraph_vector_int_clear(res);

  if (weights == NULL || igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Invalid or missing weight vector", IGRAPH_EINVAL);
  }

  // Pair edges with their weights
  typedef struct {
    igraph_integer_t edge_id;
    igraph_real_t weight;
  } edge_weight_t;

  edge_weight_t *edge_weight_pairs =
      (edge_weight_t *)malloc(no_of_edges * sizeof(edge_weight_t));
  if (!edge_weight_pairs) {
    IGRAPH_ERROR("Cannot allocate memory for edge-weight pairs.",
                 IGRAPH_ENOMEM);
  }

  for (igraph_integer_t i = 0; i < no_of_edges; i++) {
    edge_weight_pairs[i].edge_id = i;
    edge_weight_pairs[i].weight = VECTOR(*weights)[i];
  }

  // Sort edges by weight (non-decreasing)
  qsort(edge_weight_pairs, no_of_edges, sizeof(edge_weight_t),
        [](const void *a, const void *b) {
          edge_weight_t *ea = (edge_weight_t *)a;
          edge_weight_t *eb = (edge_weight_t *)b;
          return (ea->weight < eb->weight)   ? -1
                 : (ea->weight > eb->weight) ? 1
                                             : 0;
        });

  // Initialize Union-Find data structure
  igraph_union_find_t uf;
  IGRAPH_CHECK(igraph_union_find_init(&uf, no_of_nodes));
  IGRAPH_FINALLY(igraph_union_find_destroy, &uf);

  // Kruskal's algorithm: add edges without creating cycles
  for (igraph_integer_t i = 0; i < no_of_edges; i++) {
    igraph_integer_t edge = edge_weight_pairs[i].edge_id;
    igraph_integer_t from = IGRAPH_FROM(graph, edge);
    igraph_integer_t to = IGRAPH_TO(graph, edge);

    if (igraph_union_find_find(&uf, from) != igraph_union_find_find(&uf, to)) {
      igraph_vector_int_push_back(res, edge);  // Add edge to result
      igraph_union_find_union(&uf, from,
                              to);  // Union the sets of the two vertices
    }

    // Stop if we already have (V-1) edges in the MST
    if (igraph_vector_int_size(res) == no_of_nodes - 1) {
      break;
    }
  }

  free(edge_weight_pairs);
  igraph_union_find_destroy(&uf);
  IGRAPH_FINALLY_CLEAN(2);  // Clean up both uf and sorted_edges

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
