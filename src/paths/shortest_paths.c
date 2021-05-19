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
static int igraph_i_average_path_length_unweighted(
        const igraph_t *graph,
        igraph_real_t *res,
        igraph_real_t *unconnected_pairs, /* if not NULL, will be set to the no. of non-connected ordered vertex pairs */
        const igraph_bool_t directed,
        const igraph_bool_t invert, /* average inverse distances instead of distances */
        const igraph_bool_t unconn  /* average over connected pairs instead of all pairs */)
{
    long int no_of_nodes = igraph_vcount(graph);
    long int source, j, n;
    long int *already_added;
    igraph_real_t no_of_pairs = no_of_nodes > 0 ? no_of_nodes * (no_of_nodes - 1.0) : 0.0; /* no. of ordered vertex pairs */
    igraph_real_t no_of_conn_pairs = 0.0; /* no. of ordered pairs between which there is a path */

    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_adjlist_t allneis;

    *res = 0;
    already_added = IGRAPH_CALLOC(no_of_nodes, long int);
    if (already_added == 0) {
        IGRAPH_ERROR("Average path length calculation failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &allneis,
        directed ? IGRAPH_OUT : IGRAPH_ALL,
        IGRAPH_LOOPS, IGRAPH_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    for (source = 0; source < no_of_nodes; source++) {
        IGRAPH_CHECK(igraph_dqueue_push(&q, source));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        already_added[source] = source + 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);

            neis = igraph_adjlist_get(&allneis, actnode);
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
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
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            }
        } /* while !igraph_dqueue_empty */
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
    igraph_dqueue_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Computes the average of pairwise distances (used for igraph_average_path_length_dijkstra),
 * or of inverse pairwise distances (used for igraph_global_efficiency), in an unweighted graph.
 * Uses Dijkstra's algorithm, therefore all weights must be non-negative.
 */
static int igraph_i_average_path_length_dijkstra(
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
         computation, as IGRAPH_FINITE() might involve a function call
         and we want to spare that. -1 will denote infinity instead.
    */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_2wheap_t Q;
    igraph_lazy_inclist_t inclist;
    long int source, j;
    igraph_real_t no_of_pairs;
    igraph_real_t no_of_conn_pairs = 0.0; /* no. of ordered pairs between which there is a path */

    if (!weights) {
        return igraph_i_average_path_length_unweighted(graph, res, unconnected_pairs, directed, invert, unconn);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%ld) does not match the number of edges (%ld).",
                      IGRAPH_EINVAL, igraph_vector_size(weights), no_of_edges);
    }
    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
        }
        else if (igraph_is_nan(min)) {
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
            long int minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            long int nlen;

            if (minnei != source) {
                if (invert) {
                    *res += 1.0/(mindist - 1.0);
                } else {
                    *res += mindist - 1.0;
                }
                no_of_conn_pairs += 1;
            }

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_lazy_inclist_get(&inclist, (igraph_integer_t) minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int tto = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curdist = active ? -igraph_2wheap_get(&Q, tto) : 0.0;
                if (!has) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
                } else if (altdist < curdist) {
                    /* This is a shorter path */
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
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
 * has fewer than two vertices, or if the graph has no edges and \c unconn is set to \c TRUE,
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
 *    \c TRUE, only those vertex pairs will be included in the calculation
 *    between which there is a path. If \c FALSE, \c IGRAPH_INFINITY is returned
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

int igraph_average_path_length(const igraph_t *graph,
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
 * has fewer than two vertices, or if the graph has no edges and \c unconn is set to \c TRUE,
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
 *       version, \ref igraph_average_path_length() is called.
 * \param directed Boolean, whether to consider directed paths.
 *    Ignored for undirected graphs.
 * \param unconn If \c TRUE, only those pairs are considered for the calculation
 *    between which there is a path. If \c FALSE, \c IGRAPH_INFINITY is returned
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

int igraph_average_path_length_dijkstra(const igraph_t *graph,
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
 *       the global efficiency.
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

int igraph_global_efficiency(const igraph_t *graph, igraph_real_t *res,
                             const igraph_vector_t *weights,
                             igraph_bool_t directed)
{
    return igraph_i_average_path_length_dijkstra(graph, res, NULL, weights, directed, /* invert= */ 1, /* unconn= */ 0);
}


/****************************/
/***** Local efficiency *****/
/****************************/

static int igraph_i_local_efficiency_unweighted(
        const igraph_t *graph,
        const igraph_adjlist_t *adjlist,
        igraph_dqueue_t *q,
        long int *already_counted,
        igraph_vector_t *vertex_neis,
        igraph_vector_char_t *nei_mask,
        igraph_real_t *res,
        igraph_integer_t vertex,
        igraph_neimode_t mode)
{

    long int no_of_nodes = igraph_vcount(graph);
    long int vertex_neis_size;
    long int neighbor_count; /* unlike 'vertex_neis_size', 'neighbor_count' does not count self-loops and multi-edges */
    long int i, j;

    igraph_dqueue_clear(q);
    memset(already_counted, 0, no_of_nodes * sizeof(long int));

    IGRAPH_CHECK(igraph_neighbors(graph, vertex_neis, vertex, mode));
    vertex_neis_size = igraph_vector_size(vertex_neis);

    igraph_vector_char_fill(nei_mask, 0);
    neighbor_count = 0;
    for (i=0; i < vertex_neis_size; ++i) {
        long int v = VECTOR(*vertex_neis)[i];
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
        long int source = VECTOR(*vertex_neis)[i];
        long int reached = 0;

        IGRAPH_ALLOW_INTERRUPTION();

        if (source == vertex)
            continue;

        if (VECTOR(*nei_mask)[source] == 2)
            continue;
        VECTOR(*nei_mask)[source] = 2; /* mark neighbour as already processed */

        IGRAPH_CHECK(igraph_dqueue_push(q, source));
        IGRAPH_CHECK(igraph_dqueue_push(q, 0));
        already_counted[source] = source + 1;

        while (!igraph_dqueue_empty(q)) {
            igraph_vector_int_t *act_neis;
            long int act_neis_size;
            long int act = (long int) igraph_dqueue_pop(q);
            long int actdist = (long int) igraph_dqueue_pop(q);

            if (act != source && VECTOR(*nei_mask)[act]) {
                *res += 1.0 / actdist;
                reached++;
                if (reached == neighbor_count) {
                    igraph_dqueue_clear(q);
                    break;
                }
            }

            act_neis      = igraph_adjlist_get(adjlist, act);
            act_neis_size = igraph_vector_int_size(act_neis);
            for (j = 0; j < act_neis_size; j++) {
                long int neighbor = (long int) VECTOR(*act_neis)[j];

                if (neighbor == vertex || already_counted[neighbor] == i + 1)
                    continue;

                already_counted[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_push(q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(q, actdist + 1));
            }
        }
    }

    *res /= neighbor_count * (neighbor_count - 1.0);

    return IGRAPH_SUCCESS;
}

static int igraph_i_local_efficiency_dijkstra(
        const igraph_t *graph,
        igraph_lazy_inclist_t *inclist,
        igraph_2wheap_t *Q,
        igraph_vector_t *vertex_neis,
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
         computation, as IGRAPH_FINITE() might involve a function call
         and we want to spare that. -1 will denote infinity instead.
    */

    long int i, j;
    long int vertex_neis_size;
    long int neighbor_count; /* unlike 'inc_edges_size', 'neighbor_count' does not count self-loops or multi-edges */

    IGRAPH_CHECK(igraph_neighbors(graph, vertex_neis, vertex, mode));
    vertex_neis_size = igraph_vector_size(vertex_neis);

    igraph_vector_char_fill(nei_mask, 0);
    neighbor_count = 0;
    for (i=0; i < vertex_neis_size; ++i) {
        long int v = VECTOR(*vertex_neis)[i];
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
        long int source = VECTOR(*vertex_neis)[i];
        long int reached = 0;

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
            long int minnei = igraph_2wheap_max_index(Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(Q);
            igraph_vector_int_t *neis;
            long int nlen;

            if (minnei != source && VECTOR(*nei_mask)[minnei]) {
                *res += 1.0/(mindist - 1.0);
                reached++;
                if (reached == neighbor_count) {
                    igraph_2wheap_clear(Q);
                    break;
                }
            }

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_lazy_inclist_get(inclist, (igraph_integer_t) minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                igraph_real_t altdist, curdist;
                igraph_bool_t active, has;
                long int edge = (long int) VECTOR(*neis)[j];
                long int tto = IGRAPH_OTHER(graph, edge, minnei);

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
                    IGRAPH_CHECK(igraph_2wheap_modify(Q, tto, -altdist));
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
 *       \ref igraph_average_path_length() is called.
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

int igraph_local_efficiency(const igraph_t *graph, igraph_vector_t *res,
                            const igraph_vs_t vids,
                            const igraph_vector_t *weights,
                            igraph_bool_t directed, igraph_neimode_t mode)
{
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int nodes_to_calc; /* no. of vertices includes in computation */
    igraph_vit_t vit;
    igraph_vector_t vertex_neis;
    igraph_vector_char_t nei_mask;
    long int i;

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
    IGRAPH_VECTOR_INIT_FINALLY(&vertex_neis, 0);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    nodes_to_calc = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));

    if (! weights) /* unweighted case */
    {
        long int *already_counted;
        igraph_adjlist_t adjlist;
        igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;

        already_counted = IGRAPH_CALLOC(no_of_nodes, long int);
        if (already_counted == 0) {
            IGRAPH_ERROR("Local efficiency calculation failed", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, already_counted);

        IGRAPH_CHECK(igraph_adjlist_init(
            graph, &adjlist,
            directed ? IGRAPH_OUT : IGRAPH_ALL,
            IGRAPH_LOOPS, IGRAPH_MULTIPLE
        ));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

        IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

        for (IGRAPH_VIT_RESET(vit), i=0;
             ! IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), i++)
        {
            IGRAPH_CHECK(igraph_i_local_efficiency_unweighted(
                             graph, &adjlist,
                             &q, already_counted, &vertex_neis, &nei_mask,
                             &(VECTOR(*res)[i]), IGRAPH_VIT_GET(vit), mode));
        }

        igraph_dqueue_destroy(&q);
        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FREE(already_counted);
        IGRAPH_FINALLY_CLEAN(3);
    }
    else /* weighted case */
    {
        igraph_lazy_inclist_t inclist;
        igraph_2wheap_t Q;

        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match the number of edges", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERROR("Weight vector must be non-negative", IGRAPH_EINVAL);
            }
            else if (igraph_is_nan(min)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
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
    igraph_vector_destroy(&vertex_neis);
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
 *    If a null pointer is given, all weights are assumed to be 1.
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

int igraph_average_local_efficiency(const igraph_t *graph, igraph_real_t *res,
                                    const igraph_vector_t *weights,
                                    igraph_bool_t directed, igraph_neimode_t mode)
{
    long int no_of_nodes = igraph_vcount(graph);
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
 * The diameter of a graph is the length of the longest shortest path it has.
 * This function computes both the diameter, as well as the corresponding path.
 * The diameter of the null graph is considered be infinity by convention.
 *
 * If the graph has no vertices, \c IGRAPH_NAN is returned.
 *
 * \param graph The graph object.
 * \param pres Pointer to a real number, if not \c NULL then it will contain
 *        the diameter (the actual distance).
 * \param pfrom Pointer to an integer, if not \c NULL it will be set to the
 *        source vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param pto Pointer to an integer, if not \c NULL it will be set to the
 *        target vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param path Pointer to an initialized vector. If not \c NULL the actual
 *        longest geodesic path will be stored here. The vector will be
 *        resized as needed.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c TRUE the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is returned.
 * \return Error code:
 *         \c IGRAPH_ENOMEM, not enough memory for
 *         temporary data.
 *
 * Time complexity: O(|V||E|), the
 * number of vertices times the number of edges.
 *
 * \sa \ref igraph_diameter_dijkstra()
 *
 * \example examples/simple/igraph_diameter.c
 */

int igraph_diameter(const igraph_t *graph, igraph_real_t *pres,
                    igraph_integer_t *pfrom, igraph_integer_t *pto,
                    igraph_vector_t *path,
                    igraph_bool_t directed, igraph_bool_t unconn) {

    long int no_of_nodes = igraph_vcount(graph);
    long int i, j, n;
    long int *already_added;
    long int nodes_reached;
    long int from = 0, to = 0;
    igraph_real_t res = 0;

    igraph_dqueue_t q = IGRAPH_DQUEUE_NULL;
    igraph_vector_int_t *neis;
    igraph_neimode_t dirmode;
    igraph_adjlist_t allneis;

    /* See https://github.com/igraph/igraph/issues/1538#issuecomment-724071857
     * for why we return NaN for the null graph. */
    if (no_of_nodes == 0) {
        if (pres) {
            *pres = IGRAPH_NAN;
        }
        if (path) {
            igraph_vector_clear(path);
        }
        if (pfrom) {
            *pfrom = -1;
        }
        if (pto) {
            *pto = -1;
        }
        return IGRAPH_SUCCESS;
    }

    if (directed) {
        dirmode = IGRAPH_OUT;
    } else {
        dirmode = IGRAPH_ALL;
    }
    already_added = IGRAPH_CALLOC(no_of_nodes, long int);
    if (already_added == 0) {
        IGRAPH_ERROR("diameter failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, dirmode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

    for (i = 0; i < no_of_nodes; i++) {
        nodes_reached = 1;
        IGRAPH_CHECK(igraph_dqueue_push(&q, i));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        already_added[i] = i + 1;

        IGRAPH_PROGRESS("Diameter: ", 100.0 * i / no_of_nodes, NULL);

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actdist = (long int) igraph_dqueue_pop(&q);
            if (actdist > res) {
                res = actdist;
                from = i;
                to = actnode;
            }

            neis = igraph_adjlist_get(&allneis, actnode);
            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                long int neighbor = (long int) VECTOR(*neis)[j];
                if (already_added[neighbor] == i + 1) {
                    continue;
                }
                already_added[neighbor] = i + 1;
                nodes_reached++;
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
            }
        } /* while !igraph_dqueue_empty */

        /* not connected, return IGRAPH_INFINITY */
        if (nodes_reached != no_of_nodes && !unconn) {
            res = IGRAPH_INFINITY;
            from = -1;
            to = -1;
            break;
        }
    } /* for i<no_of_nodes */

    IGRAPH_PROGRESS("Diameter: ", 100.0, NULL);

    /* return the requested info */
    if (pres != 0) {
        *pres = res;
    }
    if (pfrom != 0) {
        *pfrom = (igraph_integer_t) from;
    }
    if (pto != 0) {
        *pto = (igraph_integer_t) to;
    }
    if (path != 0) {
        if (! igraph_finite(res)) {
            igraph_vector_clear(path);
        } else {
            igraph_vector_ptr_t tmpptr;
            igraph_vector_ptr_init(&tmpptr, 1);
            IGRAPH_FINALLY(igraph_vector_ptr_destroy, &tmpptr);
            VECTOR(tmpptr)[0] = path;
            IGRAPH_CHECK(igraph_get_shortest_paths(graph, &tmpptr, 0,
                                                   (igraph_integer_t) from,
                                                   igraph_vss_1((igraph_integer_t)to),
                                                   dirmode, 0, 0));
            igraph_vector_ptr_destroy(&tmpptr);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /* clean */
    IGRAPH_FREE(already_added);
    igraph_dqueue_destroy(&q);
    igraph_adjlist_destroy(&allneis);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_diameter_dijkstra
 * \brief Calculates the weighted diameter of a graph using Dijkstra's algorithm.
 *
 * This function computes the weighted diameter of a graph.
 *
 * If the graph has no vertices, \c IGRAPH_NAN is returned.
 *
 * \param graph The input graph, can be directed or undirected.
 * \param pres Pointer to a real number, if not \c NULL then it will contain
 *        the diameter (the actual distance).
 * \param pfrom Pointer to an integer, if not \c NULL it will be set to the
 *        source vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param pto Pointer to an integer, if not \c NULL it will be set to the
 *        target vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param path Pointer to an initialized vector. If not \c NULL the actual
 *        longest geodesic path will be stored here. The vector will be
 *        resized as needed.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c TRUE the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is
 *        returned.
 * \return Error code.
 *
 * Time complexity: O(|V||E|*log|E|), |V| is the number of vertices,
 * |E| is the number of edges.
 *
 * \sa \ref igraph_diameter()
 */


int igraph_diameter_dijkstra(const igraph_t *graph,
                             const igraph_vector_t *weights,
                             igraph_real_t *pres,
                             igraph_integer_t *pfrom,
                             igraph_integer_t *pto,
                             igraph_vector_t *path,
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
       - we don't use IGRAPH_INFINITY during the computation, as IGRAPH_FINITE()
         might involve a function call and we want to spare that. -1 will denote
         infinity instead.
    */

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);

    igraph_2wheap_t Q;
    igraph_inclist_t inclist;
    long int source, j;
    igraph_neimode_t dirmode = directed ? IGRAPH_OUT : IGRAPH_ALL;

    long int from = -1, to = -1;
    igraph_real_t res = 0;
    long int nodes_reached = 0;

    /* See https://github.com/igraph/igraph/issues/1538#issuecomment-724071857
     * for why we return NaN for the null graph. */
    if (no_of_nodes == 0) {
        if (pres) {
            *pres = IGRAPH_NAN;
        }
        if (path) {
            igraph_vector_clear(path);
        }
        if (pfrom) {
            *pfrom = -1;
        }
        if (pto) {
            *pto = -1;
        }
        return IGRAPH_SUCCESS;
    }

    if (!weights) {
        igraph_real_t diameter;
        IGRAPH_CHECK(igraph_diameter(graph, &diameter, pfrom, pto, path, directed, unconn));
        if (pres) {
            *pres = diameter;
        }
        return IGRAPH_SUCCESS;
    }

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERROR("Weight vector must be non-negative", IGRAPH_EINVAL);
        }
        else if (igraph_is_nan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
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
            long int minnei = igraph_2wheap_max_index(&Q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
            igraph_vector_int_t *neis;
            long int nlen;

            if (mindist > res) {
                res = mindist; from = source; to = minnei;
            }
            nodes_reached++;

            /* Now check all neighbors of 'minnei' for a shorter path */
            neis = igraph_inclist_get(&inclist, minnei);
            nlen = igraph_vector_int_size(neis);
            for (j = 0; j < nlen; j++) {
                long int edge = (long int) VECTOR(*neis)[j];
                long int tto = IGRAPH_OTHER(graph, edge, minnei);
                igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
                igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
                igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
                igraph_real_t curdist = active ? -igraph_2wheap_get(&Q, tto) : 0.0;

                if (!has) {
                    /* First finite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
                } else if (altdist < curdist) {
                    /* A shorter path */
                    IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
                }
            }

        } /* !igraph_2wheap_empty(&Q) */

        /* not connected, return infinity */
        if (nodes_reached != no_of_nodes && !unconn) {
            res = IGRAPH_INFINITY;
            from = to = -1;
            break;
        }

    } /* source < no_of_nodes */

    /* Compensate for the +1 that we have added to distances */
    res -= 1;

    igraph_inclist_destroy(&inclist);
    igraph_2wheap_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_PROGRESS("Weighted diameter: ", 100.0, NULL);

    if (pres) {
        *pres = res;
    }
    if (pfrom) {
        *pfrom = (igraph_integer_t) from;
    }
    if (pto) {
        *pto = (igraph_integer_t) to;
    }
    if (path) {
        if (!igraph_finite(res)) {
            igraph_vector_clear(path);
        } else {
            igraph_vector_ptr_t tmpptr;
            igraph_vector_ptr_init(&tmpptr, 1);
            IGRAPH_FINALLY(igraph_vector_ptr_destroy, &tmpptr);
            VECTOR(tmpptr)[0] = path;
            IGRAPH_CHECK(igraph_get_shortest_paths_dijkstra(graph,
                         /*vertices=*/ &tmpptr, /*edges=*/ 0,
                         (igraph_integer_t) from,
                         igraph_vss_1((igraph_integer_t) to),
                         weights, dirmode, /*predecessors=*/ 0,
                         /*inbound_edges=*/ 0));
            igraph_vector_ptr_destroy(&tmpptr);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return 0;
}
