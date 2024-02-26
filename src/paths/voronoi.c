/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_nongraph.h"
#include "igraph_random.h"

#include "core/indheap.h"
#include "core/interruption.h"

static igraph_error_t igraph_i_voronoi(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_vector_t *mindist,
        const igraph_vector_int_t *generators,
        igraph_neimode_t mode,
        igraph_voronoi_tiebreaker_t tiebreaker) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_generators = igraph_vector_int_size(generators);
    igraph_adjlist_t al;
    igraph_dqueue_int_t q;

    igraph_vector_int_t already_counted;

    /* tie_count[vid] is the number of generators that vid is an equal distance from.
     * This value is needed to pick one of these generators uniformly at random
     * without needing to store all of them. */
    igraph_vector_int_t tie_count;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&already_counted, no_of_nodes);

    if (tiebreaker == IGRAPH_VORONOI_RANDOM) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&tie_count, no_of_nodes);
        RNG_BEGIN();
    }

    IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
    igraph_vector_int_fill(membership, -1);

    IGRAPH_CHECK(igraph_vector_resize(mindist, no_of_nodes));
    igraph_vector_fill(mindist, IGRAPH_INFINITY);

    /* Loop through all generator points and compute shortest paths to all other vertices.
     * As we go, we keep track of the shortest distance from any generator to each vertex
     * in 'mindist'. If we find that the distance from the current generator to a vertex
     * is shorter than what was recorded so far in 'mindist', we update 'mindist' and
     * assign that vertex to the current generator.
     */
    for (igraph_integer_t i=0; i < no_of_generators; i++) {
        igraph_integer_t g = VECTOR(*generators)[i];

        IGRAPH_ALLOW_INTERRUPTION();

        /* BFS-based unweighted shortest path implementation */

        igraph_dqueue_int_clear(&q);

        VECTOR(already_counted)[g] = i+1;
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, g));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t vid  = igraph_dqueue_int_pop(&q);
            igraph_integer_t dist = igraph_dqueue_int_pop(&q);

            /* Attention! This must be igraph_real_t, not igraph_integer_t
             * because later it will be compared with another igraph_real_t
             * whose value may be infinite. */
            igraph_real_t md = VECTOR(*mindist)[vid];

            if (dist > md) {
                /* This vertex is reachable at a shorter distance from
                 * another generator. Thus all its descendants in the shortest
                 * path tree are also reachable at a shorter distance from the
                 * other generator than from the current one. Therefore
                 * we do not need to search further from here. */
                continue;
            } else if (dist < md) {
                /* This vertex is closest to the current generator so far.
                 * Assign it to the current partition. */
                VECTOR(*mindist)[vid] = dist;
                VECTOR(*membership)[vid] = i;
                if (tiebreaker == IGRAPH_VORONOI_RANDOM) {
                    VECTOR(tie_count)[vid] = 1;
                }
            } else { /* md == dist, we have a tie */
                switch (tiebreaker) {
                case IGRAPH_VORONOI_FIRST:
                    /* Never replace existing generator assignment. */
                    break;
                case IGRAPH_VORONOI_LAST:
                    /* Always replace existing generator assignment. */
                    VECTOR(*membership)[vid] = i;
                    break;
                case IGRAPH_VORONOI_RANDOM:
                    /* We replace the membership assignment with probability 1/k upon
                     * encountering the kth same-distance generator. This ensures
                     * that one of these generators is selected uniformly at random. */
                    VECTOR(tie_count)[vid]++;
                    if (RNG_UNIF01() < 1.0 / VECTOR(tie_count)[vid]) {
                        VECTOR(*membership)[vid] = i;
                    }
                    break;
                }
            }

            igraph_vector_int_t *neis = igraph_adjlist_get(&al, vid);
            igraph_integer_t nei_count = igraph_vector_int_size(neis);
            for (igraph_integer_t j = 0; j < nei_count; j++) {
                igraph_integer_t neighbor = VECTOR(*neis)[j];
                if (VECTOR(already_counted)[neighbor] == i + 1) {
                    continue;
                }
                VECTOR(already_counted)[neighbor] = i + 1;
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, dist + 1));
            }
        }
    }

    if (tiebreaker == IGRAPH_VORONOI_RANDOM) {
        RNG_END();
        igraph_vector_int_destroy(&tie_count);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&already_counted);
    igraph_dqueue_int_destroy(&q);
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_i_voronoi_dijkstra(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_vector_t *mindist,
        const igraph_vector_int_t *generators,
        const igraph_vector_t *weights,
        igraph_neimode_t mode,
        igraph_voronoi_tiebreaker_t tiebreaker) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_generators = igraph_vector_int_size(generators);
    igraph_inclist_t il;
    igraph_2wheap_t q;

    /* tie_count[vid] is the number of generators that vid is an equal distance from.
     * We use this value to be able to randomly select one of the generators. */
    igraph_vector_int_t tie_count;

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
        } else if (isnan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &il, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &il);

    IGRAPH_CHECK(igraph_2wheap_init(&q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &q);

    if (tiebreaker == IGRAPH_VORONOI_RANDOM) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&tie_count, no_of_nodes);
        RNG_BEGIN();
    }

    IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
    igraph_vector_int_fill(membership, -1);

    IGRAPH_CHECK(igraph_vector_resize(mindist, no_of_nodes));
    igraph_vector_fill(mindist, IGRAPH_INFINITY);

    /* Loop through all generator points and compute shortest paths to all other vertices.
     * As we go, we keep track of the shortest distance from any generator to each vertex
     * in 'mindist'. If we find that the distance from the current generator to a vertex
     * is shorter than what was recorded so far in 'mindist', we update 'mindist' and
     * assign that vertex to the current generator.
     */
    for (igraph_integer_t i=0; i < no_of_generators; i++) {
        igraph_integer_t g = VECTOR(*generators)[i];

        /* Weighted shortest path implementation using Dijkstra's algorithm */

        IGRAPH_ALLOW_INTERRUPTION();

        igraph_2wheap_clear(&q);

        /* We store negative distances in the maximum heap. Since some systems
         * distinguish between -0.0 and +0.0, we need -0.0 to ensure +0.0 as
         * the final result. */
        IGRAPH_CHECK(igraph_2wheap_push_with_index(&q, g, -0.0));

        while (!igraph_2wheap_empty(&q)) {
            igraph_integer_t vid = igraph_2wheap_max_index(&q);
            igraph_real_t dist = -igraph_2wheap_deactivate_max(&q);

            igraph_real_t md = VECTOR(*mindist)[vid];

            int cmp_result = igraph_cmp_epsilon(dist, md, IGRAPH_SHORTEST_PATH_EPSILON);
            if (cmp_result > 0) { /* dist > md */
                /* This vertex is reachable at a shorter distance from
                 * another generator. Thus all its descendants in the shortest
                 * path tree are also reachable at a shorter distance from the
                 * other generator than from the current one. Therefore
                 * we do not need to search further from here. */
                continue;
            } else if (cmp_result < 0) { /* dist < md */
                /* This vertex is closest to the current generator so far.
                 * Assign it to the current partition. */
                VECTOR(*mindist)[vid] = dist;
                VECTOR(*membership)[vid] = i;
                if (tiebreaker == IGRAPH_VORONOI_RANDOM) {
                    VECTOR(tie_count)[vid] = 1;
                }
            } else { /* md == dist, we have a tie */
                switch (tiebreaker) {
                case IGRAPH_VORONOI_FIRST:
                    /* Never replace existing generator assignment. */
                    break;
                case IGRAPH_VORONOI_LAST:
                    /* Always replace existing generator assignment. */
                    VECTOR(*membership)[vid] = i;
                    break;
                case IGRAPH_VORONOI_RANDOM:
                    /* We replace the membership assignment with probability 1/k upon
                     * encountering the kth same-distance generator. This ensures
                     * that one of these generators is selected uniformly at random. */
                    VECTOR(tie_count)[vid]++;
                    if (RNG_UNIF01() < 1.0 / VECTOR(tie_count)[vid]) {
                        VECTOR(*membership)[vid] = i;
                    }
                    break;
                }
            }

            igraph_vector_int_t *inc_edges = igraph_inclist_get(&il, vid);
            igraph_integer_t inc_count = igraph_vector_int_size(inc_edges);
            for (igraph_integer_t j=0; j < inc_count; j++) {
                igraph_integer_t edge = VECTOR(*inc_edges)[j];
                igraph_real_t weight = VECTOR(*weights)[edge];

                /* Optimization: do not follow infinite-weight edges. */
                if (weight == IGRAPH_INFINITY) {
                    continue;
                }

                igraph_integer_t to = IGRAPH_OTHER(graph, edge, vid);
                igraph_real_t altdist = dist + weight;

                if (! igraph_2wheap_has_elem(&q, to)) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&q, to, -altdist));
                } else if (igraph_2wheap_has_active(&q, to)) {
                    igraph_real_t curdist = -igraph_2wheap_get(&q, to);
                    if (altdist < curdist) {
                        /* This is a shorter path */
                        igraph_2wheap_modify(&q, to, -altdist);
                    }
                }
            }
        } /* !igraph_2wheap_empty(&q) */
    }

    if (tiebreaker == IGRAPH_VORONOI_RANDOM) {
        RNG_END();
        igraph_vector_int_destroy(&tie_count);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_2wheap_destroy(&q);
    igraph_inclist_destroy(&il);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_voronoi
 * \brief Voronoi partitioning of a graph.
 *
 * \experimental
 *
 * To obtain a Voronoi partitioning of a graph, we start with a set of generator
 * vertices, which will define the partitions. Each vertex is assigned to the generator
 * vertex from (or to) which it is closest.
 *
 * </para><para>
 * This function uses a BFS search for unweighted graphs and Dijkstra's algorithm
 * for weights ones.
 *
 * \param graph The graph to partition.
 * \param membership If not \c NULL, the Voronoi partition of each vertex
 *    will be stored here. <code>membership[v]</code> will be set to the index
 *    in \p generators of the generator vertex that \c v belongs to. For vertices
 *    that are not reachable from any generator, <code>-1</code> is returned.
 * \param distances If not \c NULL, the distance of each vertex to its respective
 *    generator will be stored here. For vertices which are not reachable from
 *    any generator, \c IGRAPH_INFINITY is returned.
 * \param generators Vertex IDs of the generator vertices.
 * \param weights The edge weights, interpreted as lengths in the shortest
 *    path calculation. All weights must be non-negative.
 * \param mode In directed graphs, whether to compute distances \em from
 *    generator vertices to other vertices (\c IGRAPH_OUT), \em to generator
 *    vertices from other vertices (\c IGRAPH_IN), or ignore edge directions
 *    entirely (\c IGRAPH_ALL).
 * \param tiebreaker Controls which generator vertex to assign a vertex to
 *    when it is at equal distance from/to multiple generator vertices.
 *    \clist
 *    \cli IGRAPH_VORONOI_FIRST assign the vertex to the first generator vertex.
 *    \cli IGRAPH_VORONOI_LAST assign the vertex to the last generator vertex.
 *    \cli IGRAPH_VORONOI_RANDOM assign the vertex to a random generator vertex.
 *    \endclist
 *    Note that \c IGRAPH_VORONOI_RANDOM does not guarantee that all partitions
 *    will be contiguous. For example, if 1 and 2 are chosen as generators for the
 *    graph <code>1-3, 2-3, 3-4</code>, then 3 and 4 are at equal distance from
 *    both generators. If 3 is assigned to 2 but 4 is assigned to 1, then the
 *    partition {1, 4} will not induce a connected subgraph.
 * \return Error code.
 *
 * Time complexity: In weighted graphs, O((log |S|) |E| (log |V|) + |V|), and in
 * unweighted graphs O((log |S|) |E| + |V|), where |S| is the number of generator
 * vertices, and |V| and |E| are the number of vertices and edges in the graph.
 *
 * \sa \ref igraph_distances(), \ref igraph_distances_dijkstra().
 */
igraph_error_t igraph_voronoi(
        const igraph_t *graph,
        igraph_vector_int_t *membership,
        igraph_vector_t *distances,
        const igraph_vector_int_t *generators,
        const igraph_vector_t *weights,
        igraph_neimode_t mode,
        igraph_voronoi_tiebreaker_t tiebreaker) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t *pmembership;
    igraph_vector_int_t imembership;
    igraph_vector_t *pdistances;
    igraph_vector_t idistances;

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (tiebreaker != IGRAPH_VORONOI_FIRST &&
        tiebreaker != IGRAPH_VORONOI_LAST &&
        tiebreaker != IGRAPH_VORONOI_RANDOM) {
        IGRAPH_ERROR("Invalid tiebreaker specification during Voronoi partitioning.", IGRAPH_EINVAL);
    }

    if (! igraph_vector_int_isininterval(generators, 0, igraph_vcount(graph)-1)) {
        IGRAPH_ERROR("Invalid vertex ID given as Voronoi generator.", IGRAPH_EINVVID);
    }

    if (membership) {
        pmembership = membership;
    } else {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&imembership, no_of_nodes);
        pmembership = &imembership;
    }

    if (distances) {
        pdistances = distances;
    } else {
        IGRAPH_VECTOR_INIT_FINALLY(&idistances, no_of_nodes);
        pdistances = &idistances;
    }

    if (weights) {
        IGRAPH_CHECK(igraph_i_voronoi_dijkstra(graph, pmembership, pdistances, generators, weights, mode, tiebreaker));
    } else {
        IGRAPH_CHECK(igraph_i_voronoi(graph, pmembership, pdistances, generators, mode, tiebreaker));
    }

    if (! distances) {
        igraph_vector_destroy(&idistances);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (! membership) {
        igraph_vector_int_destroy(&imembership);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
