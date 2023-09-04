/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2020  The igraph development team <igraph@igraph.org>

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

#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_dqueue.h"
#include "igraph_memory.h"
#include "igraph_progress.h"

#include "core/indheap.h"
#include "core/interruption.h"

#include <string.h>

/*****************************************************/
/***** Average path length and global efficiency *****/
/*****************************************************/

/* Computes the average of pairwise distances (used for igraph_average_path_length),
 * or of inverse pairwise distances (used for igraph_global_efficiency), in an unweighted graph. */
static igraph_error_t igraph_i_average_path_length_unweighted(
        const igraph_t *graph,
        igraph_real_t *res,
        igraph_real_t *unconnected_pairs, /* if not NULL, will be set to the no. of non-connected ordered vertex pairs */
        const igraph_bool_t directed,
        const igraph_bool_t invert, /* average inverse distances instead of distances */
        const igraph_bool_t unconn  /* average over connected pairs instead of all pairs */)
{
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t source, j, n;
    igraph_integer_t *already_added;
    igraph_real_t no_of_pairs = no_of_nodes > 0 ? no_of_nodes * (no_of_nodes - 1.0) : 0.0; /* no. of ordered vertex pairs */
    igraph_real_t no_of_conn_pairs = 0.0; /* no. of ordered pairs between which there is a path */

    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_adjlist_t allneis;

    *res = 0;
    already_added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(already_added, "Insufficient memory for average path length.");
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &allneis,
        directed ? IGRAPH_OUT : IGRAPH_ALL,
        IGRAPH_LOOPS, IGRAPH_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    for (source = 0; source < no_of_nodes; source++) {
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, source));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        already_added[source] = source + 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);

            neis = igraph_adjlist_get(&allneis, actnode);
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t neighbor = VECTOR(*neis)[j];
                if (already_added[neighbor] == source + 1) {
                    continue;
                }
                already_added[neighbor] = source + 1;
                if (invert) {
                    *res += 1.0/(actdist + 1.0);
                } else {
                    *res += actdist + 1.0;
                }
                no_of_conn_pairs += 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
            }
        } /* while !igraph_dqueue_int_empty */
    } /* for source < no_of_nodes */


    if (no_of_pairs == 0) {
        *res = IGRAPH_NAN; /* can't average zero items */
    } else {
        if (unconn) { /* average over connected pairs */
            if (no_of_conn_pairs == 0) {
                *res = IGRAPH_NAN; /* can't average zero items */
            } else {
                *res /= no_of_conn_pairs;
            }
        } else { /* average over all pairs */
            /* no_of_conn_pairs < no_of_pairs implies that the graph is disconnected */
            if (no_of_conn_pairs < no_of_pairs && ! invert) {
                /* When invert=false, assume the distance between non-connected pairs to be infinity */
                *res = IGRAPH_INFINITY;
            } else {
                /* When invert=true, assume the inverse distance between non-connected pairs
                 * to be zero. Therefore, no special treatment is needed for disconnected graphs. */
                *res /= no_of_pairs;
            }
        }
    }

    if (unconnected_pairs)
        *unconnected_pairs = no_of_pairs - no_of_conn_pairs;

    /* clean */
    IGRAPH_FREE(already_added);
    igraph_dqueue_int_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Computes the average of pairwise distances (used for igraph_average_path_length_dijkstra),
 * or of inverse pairwise distances (used for igraph_global_efficiency), in an unweighted graph.
 * Uses Dijkstra's algorithm, therefore all weights must be non-negative.
 */
static igraph_error_t igraph_i_average_path_length_dijkstra(
        const igraph_t *graph,
        igraph_real_t *res,
        igraph_real_t *unconnected_pairs,
        const igraph_vector_t *weights,
        const igraph_bool_t directed,
        const igraph_bool_t invert, /* average inverse distances instead of distances */
        const igraph_bool_t unconn  /* average over connected pairs instead of all pairs */)
{

    /* Implementation details. This is the basic Dijkstra algorithm,
       with a binary heap. The heap is indexed, i.e. it stores not only
       the distances, but also which vertex they belong to.

       From now on we use a 2-way heap, so the distances can be queried
       directly from the heap.

       Dirty tricks:
       - the opposite of the distance is stored in the heap, as it is a
         maximum heap and we need a minimum heap.
       - we don't use IGRAPH_INFINITY in the res matrix during the
         computation, as isfinite() might involve a function call
         and we want to spare that. -1 will denote infinity instead.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    igraph_integer_t source, j;
    igraph_real_t no_of_pairs;
    igraph_real_t no_of_conn_pairs = 0.0; /* no. of ordered pairs between which there is a path */

    if (!weights) {
        return igraph_i_average_path_length_unweighted(graph, res, unconnected_pairs, directed, invert, unconn);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match the number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(weights), no_of_edges);
    }
    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
        }
        else if (isnan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    /* Avoid returning a negative zero, which would be printed as -0 in tests. */
    if (no_of_nodes > 0) {
        no_of_pairs = no_of_nodes * (no_of_nodes - 1.0);
    } else {
        no_of_pairs = 0;
    }

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_lazy_inclist_init(
        graph, &inclist, directed ? IGRAPH_OUT : IGRAPH_ALL, IGRAPH_LOOPS
    ));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    *res = 0.0;

    for (source = 0; source < no_of_nodes; ++source) {

        IGRAPH_ALLOW_INTERRUPTION();

        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, -1.0);

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            igraph_integer_t nlen;

            if (minnei != source) {
                if (invert) {
                    *res += 1.0/(mindist - 1.0);
                } else {
                    *res += mindist - 1.0;
                }
                no_of_conn_pairs += 1;
            }

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_lazy_inclist_get(&inclist, minnei);
            IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_integer_t edge = VECTOR(*neis)[j];
                igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curdist = active ? -igraph_2wheap_get(&Q, tto) : 0.0;
                if (altdist == IGRAPH_INFINITY) {
                    /* Ignore edges with positive infinite weight */
                } else if (!has) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
                } else if (altdist < curdist) {
                    /* This is a shorter path */
                    igraph_2wheap_modify(&Q, tto, -altdist);
                }
            }
        } /* !igraph_2wheap_empty(&Q) */
    } /* for source < no_of_nodes */

    if (no_of_pairs == 0) {
        *res = IGRAPH_NAN; /* can't average zero items */
    } else {
        if (unconn) { /* average over connected pairs */
            if (no_of_conn_pairs == 0) {
                *res = IGRAPH_NAN; /* can't average zero items */
            } else {
                *res /= no_of_conn_pairs;
            }
        } else { /* average over all pairs */
            /* no_of_conn_pairs < no_of_pairs implies that the graph is disconnected */
            if (no_of_conn_pairs < no_of_pairs && ! invert) {
                *res = IGRAPH_INFINITY;
            } else {
                *res /= no_of_pairs;
            }
        }
    }

    if (unconnected_pairs)
        *unconnected_pairs = no_of_pairs - no_of_conn_pairs;

    igraph_lazy_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_average_path_length
 * \brief Calculates the average unweighted shortest path length between all vertex pairs.
 *
 * </para><para>
 * If no vertex pairs can be included in the calculation, for example because the graph
 * has fewer than two vertices, or if the graph has no edges and \c unconn is set to \c true,
 * NaN is returned.
 *
 * \param graph The graph object.
 * \param res Pointer to a real number, this will contain the result.
 * \param unconn_pairs Pointer to a real number. If not a null pointer, the number of
 *    ordered vertex pairs where the second vertex is unreachable from the first one
 *    will be stored here.
 * \param directed Boolean, whether to consider directed
 *    paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *    \c true, only those vertex pairs will be included in the calculation
 *    between which there is a path. If \c false, \c IGRAPH_INFINITY is returned
 *    for disconnected graphs.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for data structures
 *
 * Time complexity: O(|V| |E|), the number of vertices times the number of edges.
 *
 * \sa \ref igraph_average_path_length_dijkstra() for the weighted version.
 *
 * \example examples/simple/igraph_average_path_length.c
 */

igraph_error_t igraph_average_path_length(const igraph_t *graph,
                               igraph_real_t *res, igraph_real_t *unconn_pairs,
                               igraph_bool_t directed, igraph_bool_t unconn)
{
    return igraph_i_average_path_length_unweighted(graph, res, unconn_pairs, directed, /* invert= */ 0, unconn);
}


/**
 * \ingroup structural
 * \function igraph_average_path_length_dijkstra
 * \brief Calculates the average weighted shortest path length between all vertex pairs.
 *
 * </para><para>
 * If no vertex pairs can be included in the calculation, for example because the graph
 * has fewer than two vertices, or if the graph has no edges and \c unconn is set to \c true,
 * NaN is returned.
 *
 * </para><para>
 * All distinct ordered vertex pairs are taken into account.
 *
 * \param graph The graph object.
 * \param res Pointer to a real number, this will contain the result.
 * \param unconn_pairs Pointer to a real number. If not a null pointer, the number of
 *    ordered vertex pairs where the second vertex is unreachable from the first one
 *    will be stored here.
 * \param weights The edge weights. All edge weights must be
 *       non-negative for Dijkstra's algorithm to work. Additionally, no
 *       edge weight may be NaN. If either case does not hold, an error
 *       is returned. If this is a null pointer, then the unweighted
 *       version, \ref igraph_average_path_length() is called. Edges with positive
 *       infinite weight are ignored.
 * \param directed Boolean, whether to consider directed paths.
 *    Ignored for undirected graphs.
 * \param unconn If \c true, only those pairs are considered for the calculation
 *    between which there is a path. If \c false, \c IGRAPH_INFINITY is returned
 *    for disconnected graphs.
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_ENOMEM
 *              not enough memory for data structures
 *         \cli IGRAPH_EINVAL
 *              invalid weight vector
 *         \endclist
 *
 * Time complexity: O(|V| |E| log|E| + |V|), where |V| is the number of
 * vertices and |E| is the number of edges.
 *
 * \sa \ref igraph_average_path_length() for a slightly faster unweighted version.
 *
 * \example examples/simple/igraph_grg_game.c
 */

igraph_error_t igraph_average_path_length_dijkstra(const igraph_t *graph,
                                        igraph_real_t *res, igraph_real_t *unconn_pairs,
                                        const igraph_vector_t *weights,
                                        igraph_bool_t directed, igraph_bool_t unconn)
{
    return igraph_i_average_path_length_dijkstra(graph, res, unconn_pairs, weights, directed, /* invert= */ 0, unconn);
}


/**
 * \ingroup structural
 * \function igraph_global_efficiency
 * \brief Calculates the global efficiency of a network.
 *
 * </para><para>
 * The global efficiency of a network is defined as the average of inverse distances
 * between all pairs of vertices: <code>E_g = 1/(N*(N-1)) sum_{i!=j} 1/d_ij</code>,
 * where N is the number of vertices.
 * The inverse distance between pairs that are not reachable from each other is considered
 * to be zero. For graphs with fewer than 2 vertices, NaN is returned.
 *
 * </para><para>
 * Reference:
 * V. Latora and M. Marchiori,
 * Efficient Behavior of Small-World Networks,
 * Phys. Rev. Lett. 87, 198701 (2001).
 * https://dx.doi.org/10.1103/PhysRevLett.87.198701
 *
 * \param graph The graph object.
 * \param res Pointer to a real number, this will contain the result.
 * \param weights The edge weights. All edge weights must be
 *       non-negative for Dijkstra's algorithm to work. Additionally, no
 *       edge weight may be NaN. If either case does not hold, an error
 *       is returned. If this is a null pointer, then the unweighted
 *       version, \ref igraph_average_path_length() is used in calculating
 *       the global efficiency. Edges with positive infinite weights are
 *       ignored.
 * \param directed Boolean, whether to consider directed paths.
 *    Ignored for undirected graphs.
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_ENOMEM
 *              not enough memory for data structures
 *         \cli IGRAPH_EINVAL
 *              invalid weight vector
 *         \endclist
 *
 * Time complexity: O(|V| |E| log|E| + |V|) for weighted graphs and
 * O(|V| |E|) for unweighted ones. |V| denotes the number of
 * vertices and |E| denotes the number of edges.
 *
 */

igraph_error_t igraph_global_efficiency(const igraph_t *graph, igraph_real_t *res,
                             const igraph_vector_t *weights,
                             igraph_bool_t directed)
{
    return igraph_i_average_path_length_dijkstra(graph, res, NULL, weights, directed, /* invert= */ 1, /* unconn= */ 0);
}


/****************************/
/***** Local efficiency *****/
/****************************/

static igraph_error_t igraph_i_local_efficiency_unweighted(
        const igraph_t *graph,
        const igraph_adjlist_t *adjlist,
        igraph_dqueue_int_t *q,
        igraph_integer_t *already_counted,
        igraph_vector_int_t *vertex_neis,
        igraph_vector_char_t *nei_mask,
        igraph_real_t *res,
        igraph_integer_t vertex,
        igraph_neimode_t mode)
{

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t vertex_neis_size;
    igraph_integer_t neighbor_count; /* unlike 'vertex_neis_size', 'neighbor_count' does not count self-loops and multi-edges */
    igraph_integer_t i, j;

    igraph_dqueue_int_clear(q);

    /* already_counted[i] is 0 iff vertex i was not reached so far, otherwise
     * it is the index of the source vertex in vertex_neis that it was reached
     * from, plus 1 */
    memset(already_counted, 0, no_of_nodes * sizeof(already_counted[0]));

    IGRAPH_CHECK(igraph_neighbors(graph, vertex_neis, vertex, mode));
    vertex_neis_size = igraph_vector_int_size(vertex_neis);

    igraph_vector_char_fill(nei_mask, 0);
    neighbor_count = 0;
    for (i=0; i < vertex_neis_size; ++i) {
        igraph_integer_t v = VECTOR(*vertex_neis)[i];
        if (v != vertex && ! VECTOR(*nei_mask)[v]) {
            VECTOR(*nei_mask)[v] = 1; /* mark as unprocessed neighbour */
            neighbor_count++;
        }
    }

    *res = 0.0;

    /* when the neighbor count is smaller than 2, we return 0.0 */
    if (neighbor_count < 2) {
        return IGRAPH_SUCCESS;
    }

    for (i=0; i < vertex_neis_size; ++i) {
        igraph_integer_t source = VECTOR(*vertex_neis)[i];
        igraph_integer_t reached = 0;

        IGRAPH_ALLOW_INTERRUPTION();

        if (source == vertex)
            continue;

        if (VECTOR(*nei_mask)[source] == 2)
            continue;

        VECTOR(*nei_mask)[source] = 2; /* mark neighbour as already processed */

        IGRAPH_CHECK(igraph_dqueue_int_push(q, source));
        IGRAPH_CHECK(igraph_dqueue_int_push(q, 0));
        already_counted[source] = i + 1;

        while (!igraph_dqueue_int_empty(q)) {
            igraph_vector_int_t *act_neis;
            igraph_integer_t act_neis_size;
            igraph_integer_t act = igraph_dqueue_int_pop(q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(q);

            if (act != source && VECTOR(*nei_mask)[act]) {
                *res += 1.0 / actdist;
                reached++;
                if (reached == neighbor_count) {
                    igraph_dqueue_int_clear(q);
                    break;
                }
            }

            act_neis      = igraph_adjlist_get(adjlist, act);
            act_neis_size = igraph_vector_int_size(act_neis);
            for (j = 0; j < act_neis_size; j++) {
                igraph_integer_t neighbor = VECTOR(*act_neis)[j];

                if (neighbor == vertex || already_counted[neighbor] == i + 1)
                    continue;

                already_counted[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(q, actdist + 1));
            }
        }
    }

    *res /= neighbor_count * (neighbor_count - 1.0);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_local_efficiency_dijkstra(
        const igraph_t *graph,
        igraph_lazy_inclist_t *inclist,
        igraph_2wheap_t *Q,
        igraph_vector_int_t *vertex_neis,
        igraph_vector_char_t *nei_mask, /* true if the corresponding node is a neighbour of 'vertex' */
        igraph_real_t *res,
        igraph_integer_t vertex,
        igraph_neimode_t mode,
        const igraph_vector_t *weights)
{

    /* Implementation details. This is the basic Dijkstra algorithm,
       with a binary heap. The heap is indexed, i.e. it stores not only
       the distances, but also which vertex they belong to.

       From now on we use a 2-way heap, so the distances can be queried
       directly from the heap.

       Dirty tricks:
       - the opposite of the distance is stored in the heap, as it is a
         maximum heap and we need a minimum heap.
       - we don't use IGRAPH_INFINITY in the res matrix during the
         computation, as isfinite() might involve a function call
         and we want to spare that. -1 will denote infinity instead.
    */

    igraph_integer_t i, j;
    igraph_integer_t vertex_neis_size;
    igraph_integer_t neighbor_count; /* unlike 'inc_edges_size', 'neighbor_count' does not count self-loops or multi-edges */

    IGRAPH_CHECK(igraph_neighbors(graph, vertex_neis, vertex, mode));
    vertex_neis_size = igraph_vector_int_size(vertex_neis);

    igraph_vector_char_fill(nei_mask, 0);
    neighbor_count = 0;
    for (i=0; i < vertex_neis_size; ++i) {
        igraph_integer_t v = VECTOR(*vertex_neis)[i];
        if (v != vertex && ! VECTOR(*nei_mask)[v]) {
            VECTOR(*nei_mask)[v] = 1; /* mark as unprocessed neighbour */
            neighbor_count++;
        }
    }

    *res = 0.0;

    /* when the neighbor count is smaller than 2, we return 0.0 */
    if (neighbor_count < 2) {
        return IGRAPH_SUCCESS;
    }

    for (i=0; i < vertex_neis_size; ++i) {
        igraph_integer_t source = VECTOR(*vertex_neis)[i];
        igraph_integer_t reached = 0;

        IGRAPH_ALLOW_INTERRUPTION();

        if (source == vertex)
            continue;

        /* avoid processing a neighbour twice in multigraphs */
        if (VECTOR(*nei_mask)[source] == 2)
            continue;
        VECTOR(*nei_mask)[source] = 2; /* mark as already processed */

        igraph_2wheap_clear(Q);
        igraph_2wheap_push_with_index(Q, source, -1.0);

        while (!igraph_2wheap_empty(Q)) {
            igraph_integer_t minnei = igraph_2wheap_max_index(Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(Q);
            igraph_vector_int_t *neis;
            igraph_integer_t nlen;

            if (minnei != source && VECTOR(*nei_mask)[minnei]) {
                *res += 1.0/(mindist - 1.0);
                reached++;
                if (reached == neighbor_count) {
                    igraph_2wheap_clear(Q);
                    break;
                }
            }

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_lazy_inclist_get(inclist, minnei);
            IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_real_t altdist, curdist;
                igraph_bool_t active, has;
                igraph_integer_t edge = VECTOR(*neis)[j];
                igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);

                if (tto == vertex)
                    continue;

                altdist = mindist + VECTOR(*weights)[edge];
                active = igraph_2wheap_has_active(Q, tto);
                has = igraph_2wheap_has_elem(Q, tto);
                curdist = active ? -igraph_2wheap_get(Q, tto) : 0.0;
                if (!has) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(Q, tto, -altdist));
                } else if (altdist < curdist) {
                    /* This is a shorter path */
                    igraph_2wheap_modify(Q, tto, -altdist);
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

    }

    *res /= neighbor_count * (neighbor_count - 1.0);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_local_efficiency
 * \brief Calculates the local efficiency around each vertex in a network.
 *
 * </para><para>
 * The local efficiency of a network around a vertex is defined as follows:
 * We remove the vertex and compute the distances (shortest path lengths) between
 * its neighbours through the rest of the network. The local efficiency around the
 * removed vertex is the average of the inverse of these distances.
 *
 * </para><para>
 * The inverse distance between two vertices which are not reachable from each other
 * is considered to be zero. The local efficiency around a vertex with fewer than two
 * neighbours is taken to be zero by convention.
 *
 * </para><para>
 * Reference:
 * I. Vragović, E. Louis, and A. Díaz-Guilera,
 * Efficiency of informational transfer in regular and complex networks,
 * Phys. Rev. E 71, 1 (2005).
 * http://dx.doi.org/10.1103/PhysRevE.71.036122
 *
 * \param graph The graph object.
 * \param res Pointer to an initialized vector, this will contain the result.
 * \param vids The vertices around which the local efficiency will be calculated.
 * \param weights The edge weights. All edge weights must be
 *       non-negative. Additionally, no edge weight may be NaN. If either
 *       case does not hold, an error is returned. If this is a null
 *       pointer, then the unweighted version,
 *       \ref igraph_average_path_length() is called. Edges with positive
 *       infinite weights are ignored.
 * \param directed Boolean, whether to consider directed paths.
 *    Ignored for undirected graphs.
 * \param mode How to determine the local neighborhood of each vertex
 *    in directed graphs. Ignored in undirected graphs.
 *         \clist
 *         \cli IGRAPH_ALL
 *              take both in- and out-neighbours;
 *              this is a reasonable default for high-level interfaces.
 *         \cli IGRAPH_OUT
 *              take only out-neighbours
 *         \cli IGRAPH_IN
 *              take only in-neighbours
 *         \endclist
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_ENOMEM
 *              not enough memory for data structures
 *         \cli IGRAPH_EINVAL
 *              invalid weight vector
 *         \endclist
 *
 * Time complexity: O(|E|^2 log|E|) for weighted graphs and
 * O(|E|^2) for unweighted ones. |E| denotes the number of edges.
 *
 * \sa \ref igraph_average_local_efficiency()
 *
 */

igraph_error_t igraph_local_efficiency(const igraph_t *graph, igraph_vector_t *res,
                            const igraph_vs_t vids,
                            const igraph_vector_t *weights,
                            igraph_bool_t directed, igraph_neimode_t mode)
{
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t nodes_to_calc; /* no. of vertices includes in computation */
    igraph_vit_t vit;
    igraph_vector_int_t vertex_neis;
    igraph_vector_char_t nei_mask;
    igraph_integer_t i;

    /* 'nei_mask' is a vector indexed by vertices. The meaning of its values is as follows:
     *   0: not a neighbour of 'vertex'
     *   1: a not-yet-processed neighbour of 'vertex'
     *   2: an already processed neighbour of 'vertex'
     *
     * Marking neighbours of already processed is necessary to avoid processing them more
     * than once in multigraphs.
     */
    IGRAPH_CHECK(igraph_vector_char_init(&nei_mask, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_char_destroy, &nei_mask);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_neis, 0);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));

    if (! weights) /* unweighted case */
    {
        igraph_integer_t *already_counted;
        igraph_adjlist_t adjlist;
        igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;

        already_counted = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
        IGRAPH_CHECK_OOM(already_counted, "Insufficient memory for local efficiency calculation.");
        IGRAPH_FINALLY(igraph_free, already_counted);

        IGRAPH_CHECK(igraph_adjlist_init(
            graph, &adjlist,
            directed ? IGRAPH_OUT : IGRAPH_ALL,
            IGRAPH_LOOPS, IGRAPH_MULTIPLE
        ));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

        for (IGRAPH_VIT_RESET(vit), i=0;
             ! IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++)
        {
            IGRAPH_CHECK(igraph_i_local_efficiency_unweighted(
                             graph, &adjlist,
                             &q, already_counted, &vertex_neis, &nei_mask,
                             &(VECTOR(*res)[i]), IGRAPH_VIT_GET(vit), mode));
        }

        igraph_dqueue_int_destroy(&q);
        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FREE(already_counted);
        IGRAPH_FINALLY_CLEAN(3);
    }
    else /* weighted case */
    {
        igraph_lazy_inclist_t inclist;
        igraph_2wheap_t Q;

        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match the number of edges.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERRORF("Weights must not be negative, got %g.", IGRAPH_EINVAL, min);
            }
            else if (isnan(min)) {
                IGRAPH_ERROR("Weights must not contain NaN values.", IGRAPH_EINVAL);
            }
        }

        IGRAPH_CHECK(igraph_lazy_inclist_init(
            graph, &inclist, directed ? IGRAPH_OUT : IGRAPH_ALL, IGRAPH_LOOPS
        ));
        IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);
        IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
        IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);

        for (IGRAPH_VIT_RESET(vit), i=0;
             ! IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++)
        {
            IGRAPH_CHECK(igraph_i_local_efficiency_dijkstra(
                             graph, &inclist,
                             &Q, &vertex_neis, &nei_mask,
                             &(VECTOR(*res)[i]), IGRAPH_VIT_GET(vit), mode, weights));
        }

        igraph_2wheap_destroy(&Q);
        igraph_lazy_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vit_destroy(&vit);
    igraph_vector_int_destroy(&vertex_neis);
    igraph_vector_char_destroy(&nei_mask);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup structural
 * \function igraph_average_local_efficiency
 * \brief Calculates the average local efficiency in a network.
 *
 * For the null graph, zero is returned by convention.
 *
 * \param graph The graph object.
 * \param res Pointer to a real number, this will contain the result.
 * \param weights The edge weights. They must be all non-negative.
 *    If a null pointer is given, all weights are assumed to be 1. Edges
 *    with positive infinite weight are ignored.
 * \param directed Boolean, whether to consider directed paths.
 *    Ignored for undirected graphs.
 * \param mode How to determine the local neighborhood of each vertex
 *    in directed graphs. Ignored in undirected graphs.
 *         \clist
 *         \cli IGRAPH_ALL
 *              take both in- and out-neighbours;
 *              this is a reasonable default for high-level interfaces.
 *         \cli IGRAPH_OUT
 *              take only out-neighbours
 *         \cli IGRAPH_IN
 *              take only in-neighbours
 *         \endclist
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_ENOMEM
 *              not enough memory for data structures
 *         \cli IGRAPH_EINVAL
 *              invalid weight vector
 *         \endclist
 *
 * Time complexity: O(|E|^2 log|E|) for weighted graphs and
 * O(|E|^2) for unweighted ones. |E| denotes the number of edges.
 *
 * \sa \ref igraph_local_efficiency()
 *
 */

igraph_error_t igraph_average_local_efficiency(const igraph_t *graph, igraph_real_t *res,
                                    const igraph_vector_t *weights,
                                    igraph_bool_t directed, igraph_neimode_t mode)
{
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_t local_eff;

    /* If there are fewer than 3 vertices, no vertex has more than one neighbour, thus all
       local efficiencies are zero. For the null graph, we return zero by convention. */
    if (no_of_nodes < 3) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&local_eff, no_of_nodes);

    IGRAPH_CHECK(igraph_local_efficiency(graph, &local_eff, igraph_vss_all(), weights, directed, mode));

    *res = igraph_vector_sum(&local_eff);
    *res /= no_of_nodes;

    igraph_vector_destroy(&local_eff);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/***************************/
/***** Graph diameter ******/
/***************************/

/**
 * \ingroup structural
 * \function igraph_diameter
 * \brief Calculates the diameter of a graph (longest geodesic).
 *
 * The diameter of a graph is the length of the longest shortest path it has,
 * i.e. the maximum eccentricity of the graph's vertices.
 * This function computes both the diameter, as well as a corresponding path.
 * The diameter of the null graph is considered be infinity by convention.
 *
 * If the graph has no vertices, \c IGRAPH_NAN is returned.
 *
 * \param graph The graph object.
 * \param res Pointer to a real number, if not \c NULL then it will contain
 *        the diameter (the actual distance).
 * \param from Pointer to an integer, if not \c NULL it will be set to the
 *        source vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param to Pointer to an integer, if not \c NULL it will be set to the
 *        target vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param vertex_path Pointer to an initialized vector. If not \c NULL the actual
 *        longest geodesic path in terms of vertices will be stored here. The vector will be
 *        resized as needed.
 * \param edge_path Pointer to an initialized vector. If not \c NULL the actual
 *        longest geodesic path in terms of edges will be stored here. The vector will be
 *        resized as needed.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c true the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is returned.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V||E|), the
 * number of vertices times the number of edges.
 *
 * \sa \ref igraph_diameter_dijkstra() for the weighted version,
 * \ref igraph_radius() for the minimum eccentricity.
 *
 * \example examples/simple/igraph_diameter.c
 */

igraph_error_t igraph_diameter(const igraph_t *graph, igraph_real_t *res,
                    igraph_integer_t *from, igraph_integer_t *to,
                    igraph_vector_int_t *vertex_path, igraph_vector_int_t *edge_path,
                    igraph_bool_t directed, igraph_bool_t unconn) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i, j, n;
    igraph_integer_t *already_added;
    igraph_integer_t nodes_reached;
    /* from/to are initialized to 0 because in a singleton graph, or in an edgeless graph
     * with unconn = true, the diameter path will be considered to consist of vertex 0 only. */
    igraph_integer_t ifrom = 0, ito = 0;
    igraph_real_t ires = 0;

    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_neimode_t dirmode;
    igraph_adjlist_t allneis;

    /* See https://github.com/igraph/igraph/issues/1538#issuecomment-724071857
     * for why we return NaN for the null graph. */
    if (no_of_nodes == 0) {
        if (res) {
            *res = IGRAPH_NAN;
        }
        if (vertex_path) {
            igraph_vector_int_clear(vertex_path);
        }
        if (edge_path) {
            igraph_vector_int_clear(edge_path);
        }
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
        return IGRAPH_SUCCESS;
    }

    if (directed) {
        dirmode = IGRAPH_OUT;
    } else {
        dirmode = IGRAPH_ALL;
    }
    already_added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(already_added, "Insufficient memory for diameter calculation.");
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, dirmode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    for (i = 0; i < no_of_nodes; i++) {
        nodes_reached = 1;
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        already_added[i] = i + 1;

        IGRAPH_PROGRESS("Diameter: ", 100.0 * i / no_of_nodes, NULL);

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            if (actdist > ires) {
                ires = actdist;
                ifrom = i;
                ito = actnode;
            }

            neis = igraph_adjlist_get(&allneis, actnode);
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t neighbor = VECTOR(*neis)[j];
                if (already_added[neighbor] == i + 1) {
                    continue;
                }
                already_added[neighbor] = i + 1;
                nodes_reached++;
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
            }
        } /* while !igraph_dqueue_int_empty */

        /* not connected, return IGRAPH_INFINITY */
        if (nodes_reached != no_of_nodes && !unconn) {
            ires = IGRAPH_INFINITY;
            ifrom = -1;
            ito = -1;
            break;
        }
    } /* for i<no_of_nodes */

    IGRAPH_PROGRESS("Diameter: ", 100.0, NULL);

    /* return the requested info */
    if (res != 0) {
        *res = ires;
    }
    if (from != 0) {
        *from = ifrom;
    }
    if (to != 0) {
        *to = ito;
    }
    if ((vertex_path) || (edge_path)) {
        if (! isfinite(ires)) {
            if (vertex_path) {
                igraph_vector_int_clear(vertex_path);
            }
            if (edge_path){
                igraph_vector_int_clear(edge_path);
            }
        } else {
            IGRAPH_CHECK(igraph_get_shortest_path(graph, vertex_path, edge_path,
                                                  ifrom, ito, dirmode));
        }
    }

    /* clean */
    IGRAPH_FREE(already_added);
    igraph_dqueue_int_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_diameter_dijkstra
 * \brief Calculates the weighted diameter of a graph using Dijkstra's algorithm.
 *
 * This function computes the weighted diameter of a graph, defined as the longest
 * weighted shortest path, or the maximum weighted eccentricity of the graph's
 * vertices. A corresponding shortest path, as well as its endpoints,
 * can also be optionally computed.
 *
 * If the graph has no vertices, \c IGRAPH_NAN is returned.
 *
 * \param graph The input graph, can be directed or undirected.
 * \param weights The edge weights of the graph. Can be \c NULL for an
 *        unweighted graph. Edges with positive infinite weight are ignored.
 * \param res Pointer to a real number, if not \c NULL then it will contain
 *        the diameter (the actual distance).
 * \param from Pointer to an integer, if not \c NULL it will be set to the
 *        source vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param to Pointer to an integer, if not \c NULL it will be set to the
 *        target vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param vertex_path Pointer to an initialized vector. If not \c NULL the actual
 *        longest geodesic path in terms of vertices will be stored here. The vector will be
 *        resized as needed.
 * \param edge_path Pointer to an initialized vector. If not \c NULL the actual
 *        longest geodesic path in terms of edges will be stored here. The vector will be
 *        resized as needed.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c true the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is
 *        returned.
 * \return Error code.
 *
 * Time complexity: O(|V||E|*log|E|), |V| is the number of vertices,
 * |E| is the number of edges.
 *
 * \sa \ref igraph_diameter() for the unweighted version,
 * \ref igraph_radius_dijkstra() for the minimum weighted eccentricity.
 */


igraph_error_t igraph_diameter_dijkstra(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             igraph_real_t *res,
                             igraph_integer_t *from,
                             igraph_integer_t *to,
                             igraph_vector_int_t *vertex_path,
                             igraph_vector_int_t *edge_path,
                             igraph_bool_t directed,
                             igraph_bool_t unconn) {

    /* Implementation details. This is the basic Dijkstra algorithm,
       with a binary heap. The heap is indexed, i.e. it stores not only
       the distances, but also which vertex they belong to.

       From now on we use a 2-way heap, so the distances can be queried
       directly from the heap.

       Dirty tricks:
       - the opposite of the distance is stored in the heap, as it is a
         maximum heap and we need a minimum heap.
       - we don't use IGRAPH_INFINITY during the computation, as isfinite()
         might involve a function call and we want to spare that. -1 will denote
         infinity instead.
    */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    igraph_2wheap_t Q;
    igraph_inclist_t inclist;
    igraph_integer_t source, j;
    igraph_neimode_t dirmode = directed ? IGRAPH_OUT : IGRAPH_ALL;

    /* from/to are initialized to 0 because in a singleton graph, or in an edgeless graph
     * with unconn = true, the diameter path will be considered to consist of vertex 0 only. */
    igraph_integer_t ifrom = 0, ito = 0;
    igraph_real_t ires = 0;
    igraph_integer_t nodes_reached = 0;

    /* See https://github.com/igraph/igraph/issues/1538#issuecomment-724071857
     * for why we return NaN for the null graph. */
    if (no_of_nodes == 0) {
        if (res) {
            *res = IGRAPH_NAN;
        }
        if (vertex_path) {
            igraph_vector_int_clear(vertex_path);
        }
        if (edge_path) {
            igraph_vector_int_clear(edge_path);
        }
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
        return IGRAPH_SUCCESS;
    }

    if (!weights) {
        igraph_real_t diameter;
        IGRAPH_CHECK(igraph_diameter(graph, &diameter, from, to, vertex_path, edge_path, directed, unconn));
        if (res) {
            *res = diameter;
        }
        return IGRAPH_SUCCESS;
    }

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") not equal to number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(weights), no_of_edges);
    }

    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
        }
        else if (isnan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, dirmode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    for (source = 0; source < no_of_nodes; source++) {

        IGRAPH_PROGRESS("Weighted diameter: ", source * 100.0 / no_of_nodes, NULL);
        IGRAPH_ALLOW_INTERRUPTION();

        igraph_2wheap_clear(&Q);
        igraph_2wheap_push_with_index(&Q, source, -1.0);

        nodes_reached = 0.0;

        while (!igraph_2wheap_empty(&Q)) {
            igraph_integer_t minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            igraph_integer_t nlen;

            if (mindist > ires) {
                ires = mindist; ifrom = source; ito = minnei;
            }
            nodes_reached++;

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_inclist_get(&inclist, minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_integer_t edge = VECTOR(*neis)[j];
                igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curdist = active ? -igraph_2wheap_get(&Q, tto) : 0.0;

                if (!has) {
                    /* First finite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
                } else if (altdist < curdist) {
                    /* A shorter path */
                    igraph_2wheap_modify(&Q, tto, -altdist);
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

        /* not connected, return infinity */
        if (nodes_reached != no_of_nodes && !unconn) {
            ires = IGRAPH_INFINITY;
            ifrom = ito = -1;
            break;
        }

    } /* source < no_of_nodes */

    /* Compensate for the +1 that we have added to distances */
    ires -= 1;

    igraph_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_PROGRESS("Weighted diameter: ", 100.0, NULL);

    if (res) {
        *res = ires;
    }
    if (from) {
        *from = ifrom;
    }
    if (to) {
        *to = ito;
    }
    if ((vertex_path) || (edge_path)) {
        if (!isfinite(ires)) {
            if (vertex_path){
                igraph_vector_int_clear(vertex_path);
            }
            if (edge_path) {
                igraph_vector_int_clear(edge_path);
            }
        } else {
            IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph,
                            /*vertices=*/ vertex_path, /*edges=*/ edge_path,
                            ifrom, ito,
                            weights, dirmode));
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * Temporarily removes all edges incident on the vertex with the given ID from
 * the graph by setting the weights of these edges to infinity.
 *
 * \param graph   the graph
 * \param weights the weights of the edges of the graph
 * \param vid     the ID of the vertex to remove
 * \param edges_removed  vector that records the IDs of the edges that were
 *        "removed" (i.e. their weights were set to infinity)
 * \param eids    temporary vector that is used to retrieve the IDs of the
 *        incident edges, to make this function free of memory allocations
 */
static igraph_error_t igraph_i_semidelete_vertex(
    const igraph_t *graph, igraph_vector_t *weights,
    igraph_integer_t vid, igraph_vector_int_t *edges_removed,
    igraph_vector_int_t *eids
) {
    igraph_integer_t j, n;

    IGRAPH_CHECK(igraph_incident(graph, eids, vid, IGRAPH_ALL));

    n = igraph_vector_int_size(eids);
    for (j = 0; j < n; j++) {
        igraph_integer_t eid = VECTOR(*eids)[j];
        IGRAPH_CHECK(igraph_vector_int_push_back(edges_removed, eid));
        VECTOR(*weights)[eid] = IGRAPH_INFINITY;
    }

    return IGRAPH_SUCCESS;
}

static igraph_bool_t igraph_i_has_edge_with_infinite_weight(
    const igraph_vector_int_t* path, const igraph_vector_t* weights
) {
    igraph_integer_t i, n;

    n = weights ? igraph_vector_int_size(path) : 0;
    for (i = 0; i < n; i++) {
        igraph_integer_t edge = VECTOR(*path)[i];
        if (!isfinite(VECTOR(*weights)[edge])) {
            return true;
        }
    }

    return false;
}

static igraph_real_t igraph_i_get_total_weight_of_path(
    igraph_vector_int_t* path, const igraph_vector_t* weights
) {
    igraph_integer_t i, n = igraph_vector_int_size(path);
    igraph_real_t result;

    if (weights) {
        result = 0;
        for (i = 0; i < n; i++) {
            igraph_integer_t edge = VECTOR(*path)[i];
            result += VECTOR(*weights)[edge];
        }
    } else {
        result = n;
    }

    return result;
}

/**
 * \function igraph_get_k_shortest_paths
 * \brief k shortest paths between two vertices.
 *
 * This function returns the \p k shortest paths between two vertices, in order of
 * increasing lengths.
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * Yen, Jin Y.:
 * An algorithm for finding shortest routes from all source nodes to a given
 * destination in general networks.
 * Quarterly of Applied Mathematics. 27 (4): 526–530. (1970)
 * https://doi.org/10.1090/qam/253822
 *
 * \param graph The graph object.
 * \param weights The edge weights of the graph. Can be \c NULL for an
 *        unweighted graph. Infinite weights will be treated as missing
 *        edges.
 * \param vertex_paths Pointer to an initialized list of integer vectors, the result
 *        will be stored here in \ref igraph_vector_int_t objects. Each vector
 *        object contains the vertex IDs along the <code>k</code>th shortest path
 *        between \p from and \p to, where \c k is the vector list index. May
 *        be \c NULL if the vertex paths are not needed.
 * \param edge_paths Pointer to an initialized list of integer vectors, the result
 *        will be stored here in \ref igraph_vector_int_t objects. Each vector
 *        object contains the edge IDs along the <code>k</code>th shortest path
 *        between \p from and \p to, where \c k is the vector list index. May be
 *        \c NULL if the edge paths are not needed.
 * \param k The number of paths.
 * \param from The ID of the vertex from which the paths are calculated.
 * \param to The ID of the vertex to which the paths are calculated.
 * \param mode The type of paths to be used for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          The outgoing paths of \p from are calculated.
 *        \cli IGRAPH_IN
 *          The incoming paths of \p from are calculated.
 *        \cli IGRAPH_ALL
 *          The directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           Not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p from or \p to is an invalid vertex id.
 *        \cli IGRAPH_EINVMODE
 *           Invalid mode argument.
 *        \cli IGRAPH_EINVAL
 *           Invalid argument.
 *        \endclist
 *
 * \sa \ref igraph_get_all_simple_paths(), \ref igraph_get_shortest_paths(),
 * \ref igraph_get_shortest_paths_dijkstra()
 *
 * Time complexity:  k |V| (|V| log|V| + |E|), where |V| is the number of vertices,
 *                  and |E| is the number of edges.
 */
igraph_error_t igraph_get_k_shortest_paths(
    const igraph_t *graph, const igraph_vector_t *weights,
    igraph_vector_int_list_t *vertex_paths,
    igraph_vector_int_list_t *edge_paths,
    igraph_integer_t k, igraph_integer_t from, igraph_integer_t to,
    igraph_neimode_t mode
) {
    igraph_vector_int_list_t paths_pot; /* potential shortest paths */
    igraph_integer_t vertex_spur;
    igraph_vector_int_t path_spur, path_root, path_total, path_shortest;
    igraph_integer_t nr_edges_root, i_path_current, i_path, edge_path_root, vertex_root_del;
    igraph_integer_t i, n;
    igraph_vector_t current_weights;
    igraph_vector_int_t edges_removed;
    igraph_integer_t nr_edges = igraph_ecount(graph);
    igraph_bool_t infinite_path, already_in_potential_paths;
    igraph_vector_int_t *path_0;
    igraph_vector_int_t eids;
    igraph_real_t path_weight, shortest_path_weight;
    igraph_integer_t edge_paths_owned = 0;

    if (!igraph_is_directed(graph) && (mode == IGRAPH_IN || mode == IGRAPH_OUT)) {
        mode = IGRAPH_ALL;
    }

    if (vertex_paths) {
        igraph_vector_int_list_clear(vertex_paths);
    }

    if (!edge_paths) {
        /* We will need our own instance */
        edge_paths = IGRAPH_CALLOC(1, igraph_vector_int_list_t);
        IGRAPH_CHECK_OOM(edge_paths, "Cannot allocate vector for storing edge paths.");
        IGRAPH_FINALLY(igraph_free, edge_paths);
        edge_paths_owned = 1;

        IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(edge_paths, 0);
        edge_paths_owned = 2;
    }

    igraph_vector_int_list_clear(edge_paths);

    if (k == 0) {
        goto cleanup;
    }

    IGRAPH_CHECK(igraph_vector_int_list_resize(edge_paths, 1));
    path_0 = igraph_vector_int_list_get_ptr(edge_paths, 0);

    IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph,
                                                   NULL,
                                                   path_0,
                                                   from,
                                                   to,
                                                   weights,
                                                   mode));

    /* Check if there's a path. */
    infinite_path = igraph_i_has_edge_with_infinite_weight(path_0, weights);
    if (infinite_path || (from != to && igraph_vector_int_size(path_0) == 0)) {
        /* No path found. */
        igraph_vector_int_list_clear(edge_paths);
        goto cleanup;
    }

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&paths_pot, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&path_spur, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&path_root, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&path_total, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges_removed, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&current_weights, nr_edges);

    /* If weights are NULL we use a uniform weight vector where each edge has
     * a weight of 1. Later on, we replace the weights of removed edges with
     * infinities. Note that we work on a copy of the weight vector so the
     * original vector remains intact.
     */
    if (weights) {
        igraph_vector_update(&current_weights, weights);
    } else {
        igraph_vector_fill(&current_weights, 1);
    }

    for (i_path_current = 1; i_path_current < k; i_path_current++) {
        igraph_vector_int_t *path_previous = igraph_vector_int_list_tail_ptr(edge_paths);
        igraph_integer_t path_previous_length = igraph_vector_int_size(path_previous);
        for (nr_edges_root = 0; nr_edges_root < path_previous_length; nr_edges_root++) {
            /* Determine spur node. */
            if (mode == IGRAPH_OUT) {
                vertex_spur = IGRAPH_FROM(graph, VECTOR(*path_previous)[nr_edges_root]);
            } else if (mode == IGRAPH_IN) {
                vertex_spur = IGRAPH_TO(graph, VECTOR(*path_previous)[nr_edges_root]);
            } else {
                igraph_integer_t eid = VECTOR(*path_previous)[nr_edges_root];
                igraph_integer_t vertex_spur_1 = IGRAPH_FROM(graph, eid);
                igraph_integer_t vertex_spur_2 = IGRAPH_TO(graph, eid);
                igraph_integer_t vertex_spur_3;
                igraph_integer_t vertex_spur_4;
                if (nr_edges_root < path_previous_length-1) {
                    igraph_integer_t eid_next = VECTOR(*path_previous)[nr_edges_root + 1];
                    vertex_spur_3 = IGRAPH_FROM(graph, eid_next);
                    vertex_spur_4 = IGRAPH_TO(graph, eid_next);
                } else {
                    vertex_spur_3 = vertex_spur_4 = to;
                }
                if (vertex_spur_1 == vertex_spur_3 || vertex_spur_1 == vertex_spur_4) {
                    vertex_spur = vertex_spur_2;
                } else {
                    vertex_spur = vertex_spur_1;
                }
            }

            /* Determine root path. */
            IGRAPH_CHECK(igraph_vector_int_resize(&path_root, nr_edges_root));
            for (i = 0; i < nr_edges_root; i++) {
                VECTOR(path_root)[i] = VECTOR(*path_previous)[i];
            }

            /* Remove edges that are part of the previous shortest paths which share the same root path. */
            for (i_path = 0; i_path < i_path_current; i_path++) {
                igraph_vector_int_t *path_check = igraph_vector_int_list_get_ptr(edge_paths, i_path);
                igraph_bool_t equal = true;
                for (i = 0; i < nr_edges_root; i++) {
                    if (VECTOR(path_root)[i] != VECTOR(*path_check)[i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges_removed, VECTOR(*path_check)[nr_edges_root]));
                    VECTOR(current_weights)[VECTOR(*path_check)[nr_edges_root]] = IGRAPH_INFINITY;
                }
            }

            /* pseudocode: for each node rootPathNode in rootPath except spurNode:
             *                 remove rootPathNode from Graph;
             */
            for (edge_path_root = 0; edge_path_root < nr_edges_root; edge_path_root++) {
                if (mode == IGRAPH_OUT) {
                    vertex_root_del = IGRAPH_FROM(graph, VECTOR(path_root)[edge_path_root]);
                } else if (mode == IGRAPH_IN) {
                    vertex_root_del = IGRAPH_TO(graph, VECTOR(path_root)[edge_path_root]);
                } else {
                    igraph_integer_t eid = VECTOR(*path_previous)[edge_path_root];
                    igraph_integer_t eid_next = VECTOR(*path_previous)[edge_path_root + 1];
                    igraph_integer_t vertex_root_del_1 = IGRAPH_FROM(graph, eid);
                    igraph_integer_t vertex_root_del_2 = IGRAPH_TO(graph, eid);
                    igraph_integer_t vertex_root_del_3 = IGRAPH_FROM(graph, eid_next);
                    igraph_integer_t vertex_root_del_4 = IGRAPH_TO(graph, eid_next);
                    if (vertex_root_del_1 == vertex_root_del_3 || vertex_root_del_1 == vertex_root_del_4) {
                        vertex_root_del = vertex_root_del_2;
                    } else {
                        vertex_root_del = vertex_root_del_1;
                    }
                }
                /* Remove vertex by setting incident edges to infinity */
                IGRAPH_CHECK(igraph_i_semidelete_vertex(
                    graph, &current_weights, vertex_root_del, &edges_removed,
                    &eids
                ));
            }

            /* Determine spur path */
            IGRAPH_CHECK(igraph_get_shortest_path_dijkstra(graph,
                                                           NULL,
                                                           &path_spur,
                                                           vertex_spur,
                                                           to,
                                                           &current_weights,
                                                           mode));
            infinite_path = igraph_i_has_edge_with_infinite_weight(&path_spur, &current_weights);

            /* Add total (root + spur) path to potential paths if it's not in there yet. */
            if (!infinite_path) {
                IGRAPH_CHECK(igraph_vector_int_update(&path_total, &path_root));
                IGRAPH_CHECK(igraph_vector_int_append(&path_total, &path_spur));

                already_in_potential_paths = false;
                n = igraph_vector_int_list_size(&paths_pot);
                for (i = 0; i < n; i++) {
                    if (igraph_vector_int_all_e(&path_total, igraph_vector_int_list_get_ptr(&paths_pot, i))) {
                        already_in_potential_paths = true;
                        break;
                    }
                }

                if (!already_in_potential_paths) {
                    IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(&paths_pot, &path_total));
                }
            }

            /* Cleanup */
            n = igraph_vector_int_size(&edges_removed);
            for (i = 0; i < n; i++) {
                VECTOR(current_weights)[VECTOR(edges_removed)[i]] =
                    weights ? VECTOR(*weights)[VECTOR(edges_removed)[i]] : 1;
            }
            igraph_vector_int_clear(&edges_removed);
        }

        /* Add shortest potential path to shortest paths */
        n = igraph_vector_int_list_size(&paths_pot);
        if (n == 0) {
            break;
        }

        shortest_path_weight = igraph_i_get_total_weight_of_path(
            igraph_vector_int_list_get_ptr(&paths_pot, 0), weights
        );
        i_path = 0;
        for (i = 1; i < n; i++) {
            path_weight = igraph_i_get_total_weight_of_path(
                igraph_vector_int_list_get_ptr(&paths_pot, i), weights
            );
            if (path_weight < shortest_path_weight) {
                i_path = i;
                shortest_path_weight = path_weight;
            }
        }

        IGRAPH_CHECK(igraph_vector_int_list_remove_fast(&paths_pot, i_path, &path_shortest));
        IGRAPH_CHECK(igraph_vector_int_list_push_back(edge_paths, &path_shortest));
    }

    igraph_vector_destroy(&current_weights);
    igraph_vector_int_destroy(&eids);
    igraph_vector_int_destroy(&edges_removed);
    igraph_vector_int_destroy(&path_total);
    igraph_vector_int_destroy(&path_root);
    igraph_vector_int_destroy(&path_spur);
    igraph_vector_int_list_destroy(&paths_pot);
    IGRAPH_FINALLY_CLEAN(7);

    if (vertex_paths) {
        igraph_integer_t no_of_edge_paths = igraph_vector_int_list_size(edge_paths);

        IGRAPH_CHECK(igraph_vector_int_list_resize(vertex_paths, no_of_edge_paths));
        for (i = 0; i < no_of_edge_paths; i++) {
            igraph_vector_int_t* edge_path = igraph_vector_int_list_get_ptr(edge_paths, i);
            igraph_vector_int_t* vertex_path = igraph_vector_int_list_get_ptr(vertex_paths, i);
            IGRAPH_CHECK(igraph_vertex_path_from_edge_path(graph, from, edge_path, vertex_path, mode));
        }
    }

cleanup:
    if (edge_paths_owned >= 2) {
        igraph_vector_int_list_destroy(edge_paths);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (edge_paths_owned >= 1) {
        igraph_free(edge_paths);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
