/*
   igraph library.
   Copyright (C) 2020-2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_bitset.h"
#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_stack.h"
#include "igraph_structural.h"
#include "igraph_vector.h"
#include "igraph_vector_list.h"

#include "core/interruption.h"

/* Move vertices in order to improve the quality of a partition.
 *
 * This function considers each vertex and greedily moves it to a neighboring
 * community that maximizes the improvement in the quality of a partition.
 * Only moves that strictly improve the quality are considered.
 *
 * The vertices are examined in a queue, and initially all vertices are put in the
 * queue in a random order. Vertices are popped from the queue when they are
 * examined, and only neighbors of vertices that are moved (which are not part of
 * the cluster the vertex was moved to) are pushed to the queue again.
 *
 * The \p membership vector is used as the starting point to move around vertices,
 * and is updated in-place.
 *
 */
static igraph_error_t leiden_fastmove_vertices(
        const igraph_t *graph,
        const igraph_inclist_t *edges_per_vertex,
        const igraph_vector_t *edge_weights,
        const igraph_vector_t *vertex_out_weights,
        const igraph_vector_t *vertex_in_weights,
        const igraph_real_t resolution,
        igraph_int_t *nb_clusters,
        igraph_vector_int_t *membership,
        igraph_bool_t *changed) {

    const igraph_int_t n = igraph_vcount(graph);
    const igraph_bool_t directed = (vertex_in_weights != NULL);
    igraph_dqueue_int_t unstable_vertices;
    igraph_real_t max_diff, diff;
    igraph_bitset_t neighbor_cluster_added, vertex_is_stable;
    igraph_vector_t cluster_out_weights, cluster_in_weights;
    igraph_vector_t edge_weights_per_cluster;
    igraph_vector_int_t neighbor_clusters;
    igraph_vector_int_t vertex_order;
    igraph_vector_int_t nb_vertices_per_cluster;
    igraph_stack_int_t empty_clusters;
    igraph_int_t c, nb_neigh_clusters;
    int iter = 0;

    /* Initialize queue of unstable vertices and whether vertex is stable. Only
     * unstable vertices are in the queue. */
    IGRAPH_BITSET_INIT_FINALLY(&vertex_is_stable, n);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&unstable_vertices, n);

    /* Shuffle vertices */
    IGRAPH_CHECK(igraph_vector_int_init_range(&vertex_order, 0, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vertex_order);
    igraph_vector_int_shuffle(&vertex_order);

    /* Add to the queue */
    for (igraph_int_t i = 0; i < n; i++) {
        IGRAPH_CHECK(igraph_dqueue_int_push(&unstable_vertices, VECTOR(vertex_order)[i]));
    }

    /* Initialize cluster weights and nb vertices */
    IGRAPH_VECTOR_INIT_FINALLY(&cluster_out_weights, n);
    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&cluster_in_weights, n);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nb_vertices_per_cluster, n);
    for (igraph_int_t i = 0; i < n; i++) {
        c = VECTOR(*membership)[i];
        VECTOR(cluster_out_weights)[c] += VECTOR(*vertex_out_weights)[i];
        if (directed) {
            VECTOR(cluster_in_weights)[c] += VECTOR(*vertex_in_weights)[i];
        }
        VECTOR(nb_vertices_per_cluster)[c] += 1;
    }

    /* Initialize empty clusters */
    IGRAPH_STACK_INT_INIT_FINALLY(&empty_clusters, n);
    for (c = 0; c < n; c++) {
        if (VECTOR(nb_vertices_per_cluster)[c] == 0) {
            IGRAPH_CHECK(igraph_stack_int_push(&empty_clusters, c));
        }
    }

    /* Initialize vectors to be used in calculating differences */
    IGRAPH_VECTOR_INIT_FINALLY(&edge_weights_per_cluster, n);

    /* Initialize neighboring cluster */
    IGRAPH_BITSET_INIT_FINALLY(&neighbor_cluster_added, n);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbor_clusters, n);

    /* Iterate while the queue is not empty */
    while (!igraph_dqueue_int_empty(&unstable_vertices)) {
        igraph_int_t v = igraph_dqueue_int_pop(&unstable_vertices);
        igraph_int_t best_cluster, current_cluster = VECTOR(*membership)[v];
        igraph_int_t degree;
        igraph_vector_int_t *edges;

        /* Remove vertex from current cluster */
        VECTOR(cluster_out_weights)[current_cluster] -= VECTOR(*vertex_out_weights)[v];
        if (directed) {
            VECTOR(cluster_in_weights)[current_cluster] -= VECTOR(*vertex_in_weights)[v];
        }
        VECTOR(nb_vertices_per_cluster)[current_cluster]--;
        if (VECTOR(nb_vertices_per_cluster)[current_cluster] == 0) {
            IGRAPH_CHECK(igraph_stack_int_push(&empty_clusters, current_cluster));
        }

        /* Find out neighboring clusters */
        c = igraph_stack_int_top(&empty_clusters);
        VECTOR(neighbor_clusters)[0] = c;
        IGRAPH_BIT_SET(neighbor_cluster_added, c);
        nb_neigh_clusters = 1;

        /* Determine the edge weight to each neighboring cluster */
        edges = igraph_inclist_get(edges_per_vertex, v);
        degree = igraph_vector_int_size(edges);
        for (igraph_int_t i = 0; i < degree; i++) {
            igraph_int_t e = VECTOR(*edges)[i];
            igraph_int_t u = IGRAPH_OTHER(graph, e, v);
            if (u != v) {
                c = VECTOR(*membership)[u];
                if (!IGRAPH_BIT_TEST(neighbor_cluster_added, c)) {
                    IGRAPH_BIT_SET(neighbor_cluster_added, c);
                    VECTOR(neighbor_clusters)[nb_neigh_clusters++] = c;
                }
                VECTOR(edge_weights_per_cluster)[c] += VECTOR(*edge_weights)[e];
            }
        }

        /* Calculate maximum diff */
        best_cluster = current_cluster;
        max_diff = VECTOR(edge_weights_per_cluster)[current_cluster];
        if (directed) {
            max_diff -=
                (VECTOR(*vertex_in_weights)[v]  * VECTOR(cluster_out_weights)[current_cluster] +
                 VECTOR(*vertex_out_weights)[v] * VECTOR(cluster_in_weights)[current_cluster]) * resolution;
        } else {
            max_diff -= VECTOR(*vertex_out_weights)[v] * VECTOR(cluster_out_weights)[current_cluster] * resolution;
        }
        for (igraph_int_t i = 0; i < nb_neigh_clusters; i++) {
            c = VECTOR(neighbor_clusters)[i];
            diff = VECTOR(edge_weights_per_cluster)[c];
            if (directed) {
                diff -= (VECTOR(*vertex_out_weights)[v] * VECTOR(cluster_in_weights)[c] +
                         VECTOR(*vertex_in_weights)[v]  * VECTOR(cluster_out_weights)[c]) * resolution;
            } else {
                diff -= VECTOR(*vertex_out_weights)[v] * VECTOR(cluster_out_weights)[c] * resolution;
            }
            /* Only consider strictly improving moves.
             * Note that this is important in considering convergence.
             */
            if (diff > max_diff) {
                best_cluster = c;
                max_diff = diff;
            }
            VECTOR(edge_weights_per_cluster)[c] = 0.0;
            IGRAPH_BIT_CLEAR(neighbor_cluster_added, c);
        }

        /* Move vertex to best cluster */
        VECTOR(cluster_out_weights)[best_cluster] += VECTOR(*vertex_out_weights)[v];
        if (directed) {
            VECTOR(cluster_in_weights)[best_cluster] += VECTOR(*vertex_in_weights)[v];
        }
        VECTOR(nb_vertices_per_cluster)[best_cluster]++;
        if (best_cluster == igraph_stack_int_top(&empty_clusters)) {
            igraph_stack_int_pop(&empty_clusters);
        }

        /* Mark vertex as stable */
        IGRAPH_BIT_SET(vertex_is_stable, v);

        /* Add stable neighbours that are not part of the new cluster to the queue */
        if (best_cluster != current_cluster) {
            *changed = true;
            VECTOR(*membership)[v] = best_cluster;

            for (igraph_int_t i = 0; i < degree; i++) {
                igraph_int_t e = VECTOR(*edges)[i];
                igraph_int_t u = IGRAPH_OTHER(graph, e, v);
                if (IGRAPH_BIT_TEST(vertex_is_stable, u) && VECTOR(*membership)[u] != best_cluster) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&unstable_vertices, u));
                    IGRAPH_BIT_CLEAR(vertex_is_stable, u);
                }
            }
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
    }

    IGRAPH_CHECK(igraph_reindex_membership(membership, NULL, nb_clusters));

    igraph_vector_int_destroy(&neighbor_clusters);
    igraph_bitset_destroy(&neighbor_cluster_added);
    igraph_vector_destroy(&edge_weights_per_cluster);
    igraph_stack_int_destroy(&empty_clusters);
    igraph_vector_int_destroy(&nb_vertices_per_cluster);
    if (directed) igraph_vector_destroy(&cluster_in_weights);
    igraph_vector_destroy(&cluster_out_weights);
    igraph_vector_int_destroy(&vertex_order);
    igraph_dqueue_int_destroy(&unstable_vertices);
    igraph_bitset_destroy(&vertex_is_stable);
    if (directed) {
        IGRAPH_FINALLY_CLEAN(10);
    } else {
        IGRAPH_FINALLY_CLEAN(9);
    }

    return IGRAPH_SUCCESS;
}

/* Clean a refined membership vector.
 *
 * This function examines all vertices in \p vertex_subset and updates
 * \p refined_membership to ensure that the clusters are numbered consecutively,
 * starting from \p nb_refined_clusters. The \p nb_refined_clusters is also
 * updated itself. If C is the initial \p nb_refined_clusters and C' the
 * resulting \p nb_refined_clusters, then vertices in \p vertex_subset are numbered
 * C, C + 1, ..., C' - 1.
 */
static igraph_error_t leiden_clean_refined_membership(
        const igraph_vector_int_t* vertex_subset,
        igraph_vector_int_t *refined_membership,
        igraph_int_t* nb_refined_clusters) {

    const igraph_int_t n = igraph_vector_int_size(vertex_subset);
    igraph_vector_int_t new_cluster;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&new_cluster, n);

    /* Clean clusters. We will store the new cluster + 1 so that cluster == 0
     * indicates that no membership was assigned yet. */
    *nb_refined_clusters += 1;
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_int_t v = VECTOR(*vertex_subset)[i];
        igraph_int_t c = VECTOR(*refined_membership)[v];
        if (VECTOR(new_cluster)[c] == 0) {
            VECTOR(new_cluster)[c] = *nb_refined_clusters;
            *nb_refined_clusters += 1;
        }
    }

    /* Assign new cluster */
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_int_t v = VECTOR(*vertex_subset)[i];
        igraph_int_t c = VECTOR(*refined_membership)[v];
        VECTOR(*refined_membership)[v] = VECTOR(new_cluster)[c] - 1;
    }
    /* We used the cluster + 1, so correct */
    *nb_refined_clusters -= 1;

    igraph_vector_int_destroy(&new_cluster);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Merge vertices for a subset of the vertices. This is used to refine a partition.
 *
 * The vertices included in \p vertex_subset are assumed to be the vertices i for which
 * membership[i] = cluster_subset.
 *
 * All vertices in \p vertex_subset are initialized to a singleton partition in \p
 * refined_membership. Only singleton clusters can be merged if they are
 * sufficiently well connected to the current subgraph induced by \p
 * vertex_subset.
 *
 * We only examine each vertex once. Instead of greedily choosing the maximum
 * possible cluster to merge with, the cluster is chosen randomly among all
 * possibilities that do not decrease the quality of the partition. The
 * probability of choosing a certain cluster is proportional to exp(diff/beta).
 * For beta to 0 this converges to selecting a cluster with the maximum
 * improvement. For beta to infinity this converges to a uniform distribution
 * among all eligible clusters.
 *
 * The \p refined_membership is updated for vertex in \p vertex_subset. The number
 * of refined clusters, \p nb_refined_clusters is used to set the actual refined
 * cluster membership and is updated after this routine. Within each cluster
 * (i.e. for a given \p vertex_subset), the refined membership is initially simply
 * set to 0, ..., n - 1 (for n vertices in \p vertex_subset). However, for each \p
 * vertex_subset the refined membership should of course be unique. Hence, after
 * merging, the refined membership starts with \p nb_refined_clusters, which is
 * also updated to ensure that the resulting \p nb_refined_clusters counts all
 * refined clusters that have already been processed. See
 * leiden_clean_refined_membership for more information about
 * this aspect.
 */
static igraph_error_t leiden_merge_vertices(
        const igraph_t *graph,
        const igraph_inclist_t *edges_per_vertex,
        const igraph_vector_t *edge_weights,
        const igraph_vector_t *vertex_out_weights,
        const igraph_vector_t *vertex_in_weights,
        const igraph_vector_int_t *vertex_subset,
        const igraph_vector_int_t *membership,
        const igraph_int_t cluster_subset,
        const igraph_real_t resolution,
        const igraph_real_t beta,
        igraph_int_t *nb_refined_clusters,
        igraph_vector_int_t *refined_membership) {

    const igraph_bool_t directed = (vertex_in_weights != NULL);
    igraph_vector_int_t vertex_order;
    igraph_bitset_t non_singleton_cluster, neighbor_cluster_added;
    igraph_real_t max_diff, total_cum_trans_diff, diff;
    igraph_real_t total_vertex_out_weight = 0.0, total_vertex_in_weight = 0.0;
    const igraph_int_t n = igraph_vector_int_size(vertex_subset);
    igraph_vector_t cluster_out_weights, cluster_in_weights;
    igraph_vector_t cum_trans_diff, edge_weights_per_cluster, external_edge_weight_per_cluster_in_subset;
    igraph_vector_int_t neighbor_clusters;
    igraph_vector_int_t *edges, nb_vertices_per_cluster;
    igraph_int_t degree, nb_neigh_clusters;

    /* Initialize cluster weights */
    IGRAPH_VECTOR_INIT_FINALLY(&cluster_out_weights, n);
    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&cluster_in_weights, n);
    }

    /* Initialize number of vertices per cluster */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nb_vertices_per_cluster, n);

    /* Initialize external edge weight per cluster in subset */
    IGRAPH_VECTOR_INIT_FINALLY(&external_edge_weight_per_cluster_in_subset, n);

    /* Initialize administration for a singleton partition */
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_int_t v = VECTOR(*vertex_subset)[i];
        VECTOR(*refined_membership)[v] = i;
        VECTOR(cluster_out_weights)[i] += VECTOR(*vertex_out_weights)[v];
        total_vertex_out_weight += VECTOR(*vertex_out_weights)[v];
        if (directed) {
            VECTOR(cluster_in_weights)[i] += VECTOR(*vertex_in_weights)[v];
            total_vertex_in_weight += VECTOR(*vertex_in_weights)[v];
        }
        VECTOR(nb_vertices_per_cluster)[i] += 1;

        /* Find out neighboring clusters */
        edges = igraph_inclist_get(edges_per_vertex, v);
        degree = igraph_vector_int_size(edges);
        for (igraph_int_t j = 0; j < degree; j++) {
            igraph_int_t e = VECTOR(*edges)[j];
            igraph_int_t u = IGRAPH_OTHER(graph, e, v);
            if (u != v && VECTOR(*membership)[u] == cluster_subset) {
                VECTOR(external_edge_weight_per_cluster_in_subset)[i] += VECTOR(*edge_weights)[e];
            }
        }
    }

    /* Shuffle vertices */
    IGRAPH_CHECK(igraph_vector_int_init_copy(&vertex_order, vertex_subset));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vertex_order);
    igraph_vector_int_shuffle(&vertex_order);

    /* Initialize non singleton clusters */
    IGRAPH_BITSET_INIT_FINALLY(&non_singleton_cluster, n);

    /* Initialize vectors to be used in calculating differences */
    IGRAPH_VECTOR_INIT_FINALLY(&edge_weights_per_cluster, n);

    /* Initialize neighboring cluster */
    IGRAPH_BITSET_INIT_FINALLY(&neighbor_cluster_added, n);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbor_clusters, n);

    /* Initialize cumulative transformed difference */
    IGRAPH_VECTOR_INIT_FINALLY(&cum_trans_diff, n);

    for (igraph_int_t i = 0; i < n; i++) {
        igraph_int_t v = VECTOR(vertex_order)[i];
        igraph_int_t chosen_cluster, best_cluster, current_cluster = VECTOR(*refined_membership)[v];
        igraph_real_t vertex_weight_prod;

        if (directed) {
            vertex_weight_prod =
                VECTOR(cluster_out_weights)[current_cluster] * (total_vertex_in_weight - VECTOR(cluster_in_weights)[current_cluster]) +
                VECTOR(cluster_in_weights)[current_cluster] * (total_vertex_out_weight - VECTOR(cluster_out_weights)[current_cluster]);
        } else {
            vertex_weight_prod = VECTOR(cluster_out_weights)[current_cluster] * (total_vertex_out_weight - VECTOR(cluster_out_weights)[current_cluster]);
        }

        if (!IGRAPH_BIT_TEST(non_singleton_cluster, current_cluster) &&
            (VECTOR(external_edge_weight_per_cluster_in_subset)[current_cluster] >=
             vertex_weight_prod * resolution)) {
            /* Remove vertex from current cluster, which is then a singleton by
             * definition. */
            VECTOR(cluster_out_weights)[current_cluster] = 0.0;
            if (directed) {
                VECTOR(cluster_in_weights)[current_cluster] = 0.0;
            }
            VECTOR(nb_vertices_per_cluster)[current_cluster] = 0;

            /* Find out neighboring clusters */
            edges = igraph_inclist_get(edges_per_vertex, v);
            degree = igraph_vector_int_size(edges);

            /* Also add current cluster to ensure it can be chosen. */
            VECTOR(neighbor_clusters)[0] = current_cluster;
            IGRAPH_BIT_SET(neighbor_cluster_added, current_cluster);
            nb_neigh_clusters = 1;
            for (igraph_int_t j = 0; j < degree; j++) {
                igraph_int_t e = VECTOR(*edges)[j];
                igraph_int_t u = IGRAPH_OTHER(graph, e, v);
                if (u != v && VECTOR(*membership)[u] == cluster_subset) {
                    igraph_int_t c = VECTOR(*refined_membership)[u];
                    if (!IGRAPH_BIT_TEST(neighbor_cluster_added, c)) {
                        IGRAPH_BIT_SET(neighbor_cluster_added, c);
                        VECTOR(neighbor_clusters)[nb_neigh_clusters++] = c;
                    }
                    VECTOR(edge_weights_per_cluster)[c] += VECTOR(*edge_weights)[e];
                }
            }

            /* Calculate diffs */
            best_cluster = current_cluster;
            max_diff = 0.0;
            total_cum_trans_diff = 0.0;
            for (igraph_int_t j = 0; j < nb_neigh_clusters; j++) {
                igraph_int_t c = VECTOR(neighbor_clusters)[j];

                if (directed) {
                    vertex_weight_prod =
                        VECTOR(cluster_out_weights)[c] * (total_vertex_in_weight - VECTOR(cluster_in_weights)[c]) +
                        VECTOR(cluster_in_weights)[c] * (total_vertex_out_weight - VECTOR(cluster_out_weights)[c]);
                } else {
                    vertex_weight_prod = VECTOR(cluster_out_weights)[c] * (total_vertex_out_weight - VECTOR(cluster_out_weights)[c]);
                }

                if (VECTOR(external_edge_weight_per_cluster_in_subset)[c] >= vertex_weight_prod * resolution) {
                    diff = VECTOR(edge_weights_per_cluster)[c];
                    if (directed) {
                        diff -= (VECTOR(*vertex_out_weights)[v] * VECTOR(cluster_in_weights)[c] +
                                 VECTOR(*vertex_in_weights)[v] * VECTOR(cluster_out_weights)[c]) * resolution;
                    } else {
                        diff -= VECTOR(*vertex_out_weights)[v] * VECTOR(cluster_out_weights)[c] * resolution;
                    }


                    if (diff > max_diff) {
                        best_cluster = c;
                        max_diff = diff;
                    }

                    /* Calculate the transformed difference for sampling */
                    if (diff >= 0) {
                        total_cum_trans_diff += exp(diff / beta);
                    }

                }

                VECTOR(cum_trans_diff)[j] = total_cum_trans_diff;
                VECTOR(edge_weights_per_cluster)[c] = 0.0;
                IGRAPH_BIT_CLEAR(neighbor_cluster_added, c);
            }

            /* Determine the neighboring cluster to which the currently selected vertex
             * will be moved.
             */
            if (total_cum_trans_diff < IGRAPH_INFINITY) {
                igraph_real_t r = RNG_UNIF(0, total_cum_trans_diff);
                igraph_int_t chosen_idx;
                igraph_vector_binsearch_slice(&cum_trans_diff, r, &chosen_idx, 0, nb_neigh_clusters);
                chosen_cluster = VECTOR(neighbor_clusters)[chosen_idx];
            } else {
                chosen_cluster = best_cluster;
            }

            /* Move vertex to randomly chosen cluster */
            VECTOR(cluster_out_weights)[chosen_cluster] += VECTOR(*vertex_out_weights)[v];
            if (directed) {
                VECTOR(cluster_in_weights)[chosen_cluster] += VECTOR(*vertex_in_weights)[v];
            }
            VECTOR(nb_vertices_per_cluster)[chosen_cluster]++;

            for (igraph_int_t j = 0; j < degree; j++) {
                igraph_int_t e = VECTOR(*edges)[j];
                igraph_int_t u = IGRAPH_OTHER(graph, e, v);
                if (VECTOR(*membership)[u] == cluster_subset) {
                    if (VECTOR(*refined_membership)[u] == chosen_cluster) {
                        VECTOR(external_edge_weight_per_cluster_in_subset)[chosen_cluster] -= VECTOR(*edge_weights)[e];
                    } else {
                        VECTOR(external_edge_weight_per_cluster_in_subset)[chosen_cluster] += VECTOR(*edge_weights)[e];
                    }
                }
            }

            /* Set cluster  */
            if (chosen_cluster != current_cluster) {
                VECTOR(*refined_membership)[v] = chosen_cluster;

                IGRAPH_BIT_SET(non_singleton_cluster, chosen_cluster);
            }
        } /* end if singleton and may be merged */
    }

    IGRAPH_CHECK(leiden_clean_refined_membership(vertex_subset, refined_membership, nb_refined_clusters));

    igraph_vector_destroy(&cum_trans_diff);
    igraph_vector_int_destroy(&neighbor_clusters);
    igraph_bitset_destroy(&neighbor_cluster_added);
    igraph_vector_destroy(&edge_weights_per_cluster);
    igraph_bitset_destroy(&non_singleton_cluster);
    igraph_vector_int_destroy(&vertex_order);
    igraph_vector_destroy(&external_edge_weight_per_cluster_in_subset);
    igraph_vector_int_destroy(&nb_vertices_per_cluster);
    if (directed) igraph_vector_destroy(&cluster_in_weights);
    igraph_vector_destroy(&cluster_out_weights);
    if (directed) {
        IGRAPH_FINALLY_CLEAN(10);
    } else {
        IGRAPH_FINALLY_CLEAN(9);
    }

    return IGRAPH_SUCCESS;
}

/* Create clusters out of a membership vector.
 *
 * It is assumed that the incoming list of integer vectors is already sized
 * appropriately (i.e. it has at least as many items as the number of clusters
 * in the membership vector), and that each item in the list of integer vectors
 * is empty.
 */
static igraph_error_t leiden_get_clusters(
        const igraph_vector_int_t *membership,
        igraph_vector_int_list_t *clusters) {

    const igraph_int_t n = igraph_vector_int_size(membership);

    for (igraph_int_t i = 0; i < n; i++) {
        /* Get cluster for vertex i */
        igraph_vector_int_t *cluster = igraph_vector_int_list_get_ptr(clusters, VECTOR(*membership)[i]);

        /* Add vertex i to cluster vector */
        IGRAPH_CHECK(igraph_vector_int_push_back(cluster, i));
    }

    return IGRAPH_SUCCESS;
}

/* Aggregate the graph based on the \p refined membership while setting the
 * membership of each aggregated vertex according to the \p membership.
 *
 * Technically speaking we have that
 * aggregated_membership[refined_membership[v]] = membership[v] for each vertex v.
 *
 * The new aggregated graph is returned in \p aggregated_graph. This graph
 * object should not yet be initialized, igraph_create() is called on it, and
 * responsibility for destroying the object lies with the calling method
 *
 * The remaining results, aggregated_edge_weights, aggregate_vertex_weights and
 * aggregated_membership are all expected to be initialized.
 *
 */
static igraph_error_t leiden_aggregate(
        const igraph_t *graph,
        const igraph_inclist_t *edges_per_vertex,
        const igraph_vector_t *edge_weights,
        const igraph_vector_t *vertex_out_weights,
        const igraph_vector_t *vertex_in_weights,
        const igraph_vector_int_t *membership,
        const igraph_vector_int_t *refined_membership,
        const igraph_int_t nb_refined_clusters,
        igraph_t *aggregated_graph,
        igraph_vector_t *aggregated_edge_weights,
        igraph_vector_t *aggregated_vertex_out_weights,
        igraph_vector_t *aggregated_vertex_in_weights,
        igraph_vector_int_t *aggregated_membership) {

    const igraph_bool_t directed = (vertex_in_weights != NULL);
    igraph_vector_int_t aggregated_edges;
    igraph_vector_t edge_weight_to_cluster;
    igraph_vector_int_list_t refined_clusters;
    igraph_vector_int_t *incident_edges;
    igraph_vector_int_t neighbor_clusters;
    igraph_bitset_t neighbor_cluster_added;
    igraph_int_t c, degree, nb_neigh_clusters;

    /* Get refined clusters */
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&refined_clusters, nb_refined_clusters);
    IGRAPH_CHECK(leiden_get_clusters(refined_membership, &refined_clusters));

    /* Initialize new edges */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&aggregated_edges, 0);

    /* We clear the aggregated edge weights, we will push each new edge weight */
    igraph_vector_clear(aggregated_edge_weights);
    /* Simply resize the aggregated vertex weights and membership, they can be set directly */
    IGRAPH_CHECK(igraph_vector_resize(aggregated_vertex_out_weights, nb_refined_clusters));
    if (directed) {
        IGRAPH_CHECK(igraph_vector_resize(aggregated_vertex_in_weights, nb_refined_clusters));
    }
    IGRAPH_CHECK(igraph_vector_int_resize(aggregated_membership, nb_refined_clusters));

    IGRAPH_VECTOR_INIT_FINALLY(&edge_weight_to_cluster, nb_refined_clusters);

    /* Initialize neighboring cluster */
    IGRAPH_BITSET_INIT_FINALLY(&neighbor_cluster_added, nb_refined_clusters);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbor_clusters, nb_refined_clusters);

    /* Check per cluster */
    for (c = 0; c < nb_refined_clusters; c++) {
        igraph_vector_int_t* refined_cluster = igraph_vector_int_list_get_ptr(&refined_clusters, c);
        igraph_int_t n_c = igraph_vector_int_size(refined_cluster);
        igraph_int_t v = -1;

        /* Calculate the total edge weight to other clusters */
        VECTOR(*aggregated_vertex_out_weights)[c] = 0.0;
        if (directed) {
            VECTOR(*aggregated_vertex_in_weights)[c] = 0.0;
        }
        nb_neigh_clusters = 0;
        for (igraph_int_t i = 0; i < n_c; i++) {
            v = VECTOR(*refined_cluster)[i];
            incident_edges = igraph_inclist_get(edges_per_vertex, v);
            degree = igraph_vector_int_size(incident_edges);

            for (igraph_int_t j = 0; j < degree; j++) {
                igraph_int_t e = VECTOR(*incident_edges)[j];
                igraph_int_t u = IGRAPH_OTHER(graph, e, v);
                igraph_int_t c2 = VECTOR(*refined_membership)[u];

                if (c2 > c) {
                    if (!IGRAPH_BIT_TEST(neighbor_cluster_added, c2)) {
                        IGRAPH_BIT_SET(neighbor_cluster_added, c2);
                        VECTOR(neighbor_clusters)[nb_neigh_clusters++] = c2;
                    }
                    VECTOR(edge_weight_to_cluster)[c2] += VECTOR(*edge_weights)[e];
                }
            }

            VECTOR(*aggregated_vertex_out_weights)[c] += VECTOR(*vertex_out_weights)[v];
            if (directed) {
                VECTOR(*aggregated_vertex_in_weights)[c] += VECTOR(*vertex_in_weights)[v];
            }
        }

        /* Add actual edges from this cluster to the other clusters */
        for (igraph_int_t i = 0; i < nb_neigh_clusters; i++) {
            igraph_int_t c2 = VECTOR(neighbor_clusters)[i];

            /* Add edge */
            IGRAPH_CHECK(igraph_vector_int_push_back(&aggregated_edges, c));
            IGRAPH_CHECK(igraph_vector_int_push_back(&aggregated_edges, c2));

            /* Add edge weight */
            IGRAPH_CHECK(igraph_vector_push_back(aggregated_edge_weights, VECTOR(edge_weight_to_cluster)[c2]));

            VECTOR(edge_weight_to_cluster)[c2] = 0.0;
            IGRAPH_BIT_CLEAR(neighbor_cluster_added, c2);
        }

        VECTOR(*aggregated_membership)[c] = VECTOR(*membership)[v];

    }

    igraph_vector_int_destroy(&neighbor_clusters);
    igraph_bitset_destroy(&neighbor_cluster_added);
    igraph_vector_destroy(&edge_weight_to_cluster);
    igraph_vector_int_list_destroy(&refined_clusters);
    IGRAPH_FINALLY_CLEAN(4);

    igraph_destroy(aggregated_graph);
    IGRAPH_CHECK(igraph_create(aggregated_graph, &aggregated_edges, nb_refined_clusters,
                               directed));

    igraph_vector_int_destroy(&aggregated_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Calculate the quality of the partition.
 *
 * The quality is defined as
 *
 * 1 / 2m sum_ij (A_ij - gamma n_i n_j) d(s_i, s_j)
 *
 * for undirected graphs and as
 *
 * 1 / m sum_ij (A_ij - gamma n^out_i n^in_j) d(s_i, s_j)
 *
 * where m is the total edge weight, A_ij is the weight of edge (i, j), gamma is
 * the so-called resolution parameter, n_i is the vertex weight of vertex i, s_i is
 * the cluster of vertex i and d(x, y) = 1 if and only if x = y and 0 otherwise.
 *
 * Note that by setting n_i = k_i the degree of vertex i and dividing gamma by 2m,
 * we effectively optimize modularity. By setting n_i = 1 we optimize the
 * Constant Potts Model.
 *
 * This can be represented as a sum over clusters as
 *
 * 1 / 2m sum_c (e_c - gamma N_c^2)
 *
 * where e_c = sum_ij A_ij d(s_i, c)d(s_j, c) is the internal edge weight
 * in cluster c (or twice this value if undirected) and
 * N_c = sum_i n_i d(s_i, c) is the sum of the vertex weights inside cluster c.
 * This is how the quality is calculated in practice.
 */
static igraph_error_t leiden_quality(
        const igraph_t *graph,
        const igraph_vector_t *edge_weights,
        const igraph_vector_t *vertex_out_weights,
        const igraph_vector_t *vertex_in_weights,
        const igraph_vector_int_t *membership,
        const igraph_int_t nb_clusters,
        const igraph_real_t resolution,
        igraph_real_t *quality) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    const igraph_bool_t directed = (vertex_in_weights != NULL);
    const igraph_real_t directed_multiplier = directed ? 1.0 : 2.0;
    igraph_vector_t cluster_out_weights, cluster_in_weights;
    igraph_real_t total_edge_weight = 0.0;

    *quality = 0.0;

    for (igraph_int_t e=0; e < ecount; e++) {
        igraph_int_t from = IGRAPH_FROM(graph, e);
        igraph_int_t to = IGRAPH_TO(graph, e);
        total_edge_weight += VECTOR(*edge_weights)[e];

        /* We add the internal edge weights. */
        if (VECTOR(*membership)[from] == VECTOR(*membership)[to]) {
            *quality += directed_multiplier * VECTOR(*edge_weights)[e];
        }
    }

    /* Initialize and compute cluster weights. */

    IGRAPH_VECTOR_INIT_FINALLY(&cluster_out_weights, vcount);
    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&cluster_in_weights, vcount);
    }

    for (igraph_int_t i = 0; i < vcount; i++) {
        igraph_int_t c = VECTOR(*membership)[i];
        VECTOR(cluster_out_weights)[c] += VECTOR(*vertex_out_weights)[i];
        if (directed) {
            VECTOR(cluster_in_weights)[c] += VECTOR(*vertex_in_weights)[i];
        }
    }

    /* We subtract gamma * N^out_c * N^in_c */

    for (igraph_int_t c = 0; c < nb_clusters; c++) {
        if (directed) {
            *quality -= resolution * VECTOR(cluster_out_weights)[c] * VECTOR(cluster_in_weights)[c];
        } else {
            *quality -= resolution * VECTOR(cluster_out_weights)[c] * VECTOR(cluster_out_weights)[c];
        }
    }

    if (directed) {
        igraph_vector_destroy(&cluster_in_weights);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&cluster_out_weights);
    IGRAPH_FINALLY_CLEAN(1);

    /* We normalise by m or 2m depending on directedness */
    *quality /= (directed_multiplier * total_edge_weight);

    return IGRAPH_SUCCESS;
}

/* This is the core of the Leiden algorithm and relies on subroutines to
 * perform the three different phases: (1) local moving of vertices, (2)
 * refinement of the partition and (3) aggregation of the network based on the
 * refined partition, using the non-refined partition to create an initial
 * partition for the aggregate network.
 */
static igraph_error_t community_leiden(
        const igraph_t *graph,
        igraph_vector_t *edge_weights,
        igraph_vector_t *vertex_out_weights,
        igraph_vector_t *vertex_in_weights,
        igraph_real_t resolution,
        igraph_real_t beta,
        igraph_vector_int_t *membership,
        igraph_int_t *nb_clusters,
        igraph_real_t *quality,
        igraph_bool_t *changed) {

    const igraph_int_t n = igraph_vcount(graph);
    const igraph_bool_t directed = (vertex_in_weights != NULL);
    igraph_int_t nb_refined_clusters;
    igraph_int_t i, c;
    igraph_t aggregated_graph, *i_graph;
    igraph_vector_t aggregated_edge_weights;
    igraph_vector_t aggregated_vertex_out_weights, aggregated_vertex_in_weights;
    igraph_vector_int_t aggregated_membership;
    igraph_vector_t *i_edge_weights;
    igraph_vector_t *i_vertex_out_weights, *i_vertex_in_weights;
    igraph_vector_int_t *i_membership;
    igraph_vector_t tmp_edge_weights, tmp_vertex_out_weights, tmp_vertex_in_weights;
    igraph_vector_int_t tmp_membership;
    igraph_vector_int_t refined_membership;
    igraph_vector_int_t aggregate_vertex;
    igraph_vector_int_list_t clusters;
    igraph_inclist_t edges_per_vertex;
    igraph_bool_t continue_clustering;
    igraph_int_t level = 0;

    /* Initialize temporary weights and membership to be used in aggregation */
    IGRAPH_VECTOR_INIT_FINALLY(&tmp_edge_weights, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp_vertex_out_weights, 0);
    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&tmp_vertex_in_weights, 0);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp_membership, 0);

    /* Initialize clusters */
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&clusters, n);

    /* Initialize aggregate vertices, which initially is identical to simply the
     * vertices in the graph. */
    IGRAPH_CHECK(igraph_vector_int_init_range(&aggregate_vertex, 0, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &aggregate_vertex);

    /* Initialize refined membership */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&refined_membership, 0);

    /* Initialize aggregated graph */
    IGRAPH_CHECK(igraph_empty(&aggregated_graph, 0, directed));
    IGRAPH_FINALLY(igraph_destroy, &aggregated_graph);

    /* Initialize aggregated edge weights */
    IGRAPH_VECTOR_INIT_FINALLY(&aggregated_edge_weights, 0);

    /* Initialize aggregated vertex weights */
    IGRAPH_VECTOR_INIT_FINALLY(&aggregated_vertex_out_weights, 0);
    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&aggregated_vertex_in_weights, 0);
    }

    /* Initialize aggregated membership */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&aggregated_membership, 0);

    /* Set actual graph, weights and membership to be used. */
    i_graph = (igraph_t*)graph;
    i_edge_weights = edge_weights;
    i_vertex_out_weights = vertex_out_weights;
    i_vertex_in_weights = directed ? vertex_in_weights : NULL;
    i_membership = membership;

    /* Clean membership: ensure that cluster indices are 0 <= c < n. */
    IGRAPH_CHECK(igraph_reindex_membership(i_membership, NULL, nb_clusters));

    /* We start out with no changes, whenever a vertex is moved, this will be set to true. */
    *changed = false;
    do {

        /* Get incidence list for fast iteration */
        IGRAPH_CHECK(igraph_inclist_init( i_graph, &edges_per_vertex, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &edges_per_vertex);

        /* Move around the vertices in order to increase the quality */
        IGRAPH_CHECK(leiden_fastmove_vertices(i_graph,
                                              &edges_per_vertex,
                                              i_edge_weights,
                                              i_vertex_out_weights, i_vertex_in_weights,
                                              resolution,
                                              nb_clusters,
                                              i_membership,
                                              changed));

        /* We only continue clustering if not all clusters are represented by a
         * single vertex yet
         */
        continue_clustering = (*nb_clusters < igraph_vcount(i_graph));

        if (continue_clustering) {
            /* Set original membership */
            if (level > 0) {
                for (i = 0; i < n; i++) {
                    igraph_int_t v_aggregate = VECTOR(aggregate_vertex)[i];
                    VECTOR(*membership)[i] = VECTOR(*i_membership)[v_aggregate];
                }
            }

            /* Get vertex sets for each cluster. */
            IGRAPH_CHECK(leiden_get_clusters(i_membership, &clusters));

            /* Ensure refined membership is correct size */
            IGRAPH_CHECK(igraph_vector_int_resize(&refined_membership, igraph_vcount(i_graph)));

            /* Refine each cluster */
            nb_refined_clusters = 0;
            for (c = 0; c < *nb_clusters; c++) {
                igraph_vector_int_t* cluster = igraph_vector_int_list_get_ptr(&clusters, c);
                IGRAPH_CHECK(leiden_merge_vertices(i_graph,
                                                   &edges_per_vertex,
                                                   i_edge_weights,
                                                   i_vertex_out_weights, i_vertex_in_weights,
                                                   cluster, i_membership, c,
                                                   resolution, beta,
                                                   &nb_refined_clusters, &refined_membership));
                /* Empty cluster */
                igraph_vector_int_clear(cluster);
            }

            /* If refinement didn't aggregate anything, we aggregate on the basis of
             * the actual clustering */
            if (nb_refined_clusters >= igraph_vcount(i_graph)) {
                IGRAPH_CHECK(igraph_vector_int_update(&refined_membership, i_membership));
                nb_refined_clusters = *nb_clusters;
            }

            /* Keep track of aggregate vertex. */
            for (i = 0; i < n; i++) {
                /* Current aggregate vertex */
                igraph_int_t v_aggregate = VECTOR(aggregate_vertex)[i];
                /* New aggregate vertex */
                VECTOR(aggregate_vertex)[i] = VECTOR(refined_membership)[v_aggregate];
            }

            IGRAPH_CHECK(leiden_aggregate(
                i_graph,
                &edges_per_vertex,
                i_edge_weights,
                i_vertex_out_weights, i_vertex_in_weights,
                i_membership, &refined_membership, nb_refined_clusters,
                &aggregated_graph,
                &tmp_edge_weights,
                &tmp_vertex_out_weights, directed ? &tmp_vertex_in_weights : NULL,
                &tmp_membership));

            /* On the lowest level, the actual graph and vertex and edge weights and
             * membership are used. On higher levels, we will use the aggregated graph
             * and associated vectors.
             */
            if (level == 0) {
                /* Set actual graph, weights and membership to be used. */
                i_graph = &aggregated_graph;
                i_edge_weights = &aggregated_edge_weights;
                i_vertex_out_weights = &aggregated_vertex_out_weights;
                if (directed) {
                    i_vertex_in_weights = &aggregated_vertex_in_weights;
                }
                i_membership = &aggregated_membership;
            }

            /* Update the aggregated administration. */
            IGRAPH_CHECK(igraph_vector_update(i_edge_weights, &tmp_edge_weights));
            IGRAPH_CHECK(igraph_vector_update(i_vertex_out_weights, &tmp_vertex_out_weights));
            if (directed) {
                IGRAPH_CHECK(igraph_vector_update(i_vertex_in_weights, &tmp_vertex_in_weights));
            }
            IGRAPH_CHECK(igraph_vector_int_update(i_membership, &tmp_membership));

            level += 1;
        }

        /* We are done iterating, so we destroy the incidence list */
        igraph_inclist_destroy(&edges_per_vertex);
        IGRAPH_FINALLY_CLEAN(1);
    } while (continue_clustering);

    /* Free aggregated graph and associated vectors */
    igraph_vector_int_destroy(&aggregated_membership);
    if (directed) igraph_vector_destroy(&aggregated_vertex_in_weights);
    igraph_vector_destroy(&aggregated_vertex_out_weights);
    igraph_vector_destroy(&aggregated_edge_weights);
    igraph_destroy(&aggregated_graph);

    /* Free remaining memory */
    igraph_vector_int_destroy(&refined_membership);
    igraph_vector_int_destroy(&aggregate_vertex);
    igraph_vector_int_list_destroy(&clusters);
    igraph_vector_int_destroy(&tmp_membership);
    igraph_vector_destroy(&tmp_vertex_out_weights);
    if (directed) igraph_vector_destroy(&tmp_vertex_in_weights);
    igraph_vector_destroy(&tmp_edge_weights);

    if (directed) {
        IGRAPH_FINALLY_CLEAN(12);
    } else {
        IGRAPH_FINALLY_CLEAN(10);
    }

    /* Calculate quality */
    if (quality) {
        IGRAPH_CHECK(leiden_quality(graph,
                                    edge_weights, vertex_out_weights, vertex_in_weights,
                                    membership,
                                    *nb_clusters, resolution,
                                    quality));
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup communities
 * \function igraph_community_leiden
 * \brief Finding community structure using the Leiden algorithm.
 *
 * This function implements the Leiden algorithm for finding community
 * structure.
 *
 * </para><para>
 * It is similar to the multilevel algorithm, often called the Louvain
 * algorithm, but it is faster and yields higher quality solutions. It can
 * optimize both modularity and the Constant Potts Model, which does not suffer
 * from the resolution-limit (see Traag, Van Dooren &amp; Nesterov).
 *
 * </para><para>
 * The Leiden algorithm consists of three phases: (1) local moving of vertices, (2)
 * refinement of the partition and (3) aggregation of the network based on the
 * refined partition, using the non-refined partition to create an initial
 * partition for the aggregate network. In the local move procedure in the
 * Leiden algorithm, only vertices whose neighborhood has changed are visited. Only
 * moves that strictly improve the quality function are made. The refinement is
 * done by restarting from a singleton partition within each cluster and
 * gradually merging the subclusters. When aggregating, a single cluster may
 * then be represented by several vertices (which are the subclusters identified in
 * the refinement).
 *
 * </para><para>
 * The Leiden algorithm provides several guarantees. The Leiden algorithm is
 * typically iterated: the output of one iteration is used as the input for the
 * next iteration. At each iteration all clusters are guaranteed to be (weakly)
 * connected and well-separated. After an iteration in which nothing has
 * changed, all vertices and some parts are guaranteed to be locally optimally
 * assigned. Note that even if a single iteration did not result in any change,
 * it is still possible that a subsequent iteration might find some
 * improvement. Each iteration explores different subsets of vertices to consider
 * for moving from one cluster to another. Finally, asymptotically, all subsets
 * of all clusters are guaranteed to be locally optimally assigned. For more
 * details, please see Traag, Waltman &amp; van Eck (2019).
 *
 * </para><para>
 * The objective function being optimized is
 *
 * </para><para>
 * <code>1 / 2m sum_ij (A_ij - γ n_i n_j) δ(s_i, s_j)</code>
 *
 * </para><para>
 * in the undirected case and
 *
 * </para><para>
 * <code>1 / m sum_ij (A_ij - γ n^out_i n^in_j) δ(s_i, s_j)</code>
 *
 * </para><para>
 * in the directed case.
 * Here \c m is the total edge weight, <code>A_ij</code> is the weight of edge
 * (i, j), \c γ is the so-called resolution parameter, <code>n_i</code>
 * is the vertex weight of vertex \c i (separate out- and in-weights are used
 * with directed graphs), <code>s_i</code> is the cluster of vertex
 * \c i and <code>δ(x, y) = 1</code> if and only if <code>x = y</code> and 0
 * otherwise.
 *
 * </para><para>
 * By setting <code>n_i = k_i</code>, the degree of vertex \c i, and
 * dividing \c γ by <code>2m</code> (by \c m in the directed case), we effectively
 * obtain an expression for modularity. Hence, the standard modularity will be
 * optimized when you supply the degrees (out- and in-degrees with directed graphs)
 * as the vertex weights and by supplying as a resolution parameter
 * <code>1/(2m)</code> (<code>1/m</code> with directed graphs).
 * Use the \ref igraph_community_leiden_simple() convenience function to
 * compute vertex weights automatically for modularity maximization.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * V. A. Traag, L. Waltman, N. J. van Eck:
 * From Louvain to Leiden: guaranteeing well-connected communities.
 * Scientific Reports, 9(1), 5233 (2019).
 * http://dx.doi.org/10.1038/s41598-019-41695-z
 *
 * </para><para>
 * V. A. Traag, P. Van Dooren, and Y. Nesterov:
 * Narrow scope for resolution-limit-free community detection.
 * Phys. Rev. E 84, 016114 (2011).
 * https://doi.org/10.1103/PhysRevE.84.016114
 *
 * \param graph The input graph.
 * \param edge_weights Numeric vector containing edge weights. If \c NULL,
 *    every edge has equal weight of 1. The weights need not be non-negative.
 * \param vertex_out_weights Numeric vector containing vertex weights, or vertex
 *    out-weights for directed graphs. If \c NULL, every vertex has equal
 *    weight of 1.
 * \param vertex_in_weights Numeric vector containing vertex in-weights for
 *    directed graphs. If set to \c NULL, in-weights are assumed to be the same
 *    as out-weights, which effectively ignores edge directions.
 *    Must be \c NULL for undirected graphs.
 * \param n_iterations Iterate the core Leiden algorithm the indicated number
 *    of times. If this is a negative number, it will continue iterating until
 *    an iteration did not change the clustering. Two iterations are often
 *    sufficient, thus 2 is a reasonable default.
 * \param beta The randomness used in the refinement step when merging. A small
 *    amount of randomness (\c beta = 0.01) typically works well.
 * \param start Start from membership vector. If this is true, the optimization
 *    will start from the provided membership vector. If this is false, the
 *    optimization will start from a singleton partition.
 * \param n_iterations Iterate the core Leiden algorithm for the indicated number
 *    of times. If this is a negative number, it will continue iterating until
 *    an iteration did not change the clustering.
 * \param membership The membership vector. This is both used as the initial
 *    membership from which optimisation starts and is updated in place. It
 *    must hence be properly initialized. When finding clusters from scratch it
 *    is typically started using a singleton clustering. This can be achieved
 *    using \ref igraph_vector_int_init_range().
 * \param nb_clusters The number of clusters contained in the final \p membership.
 *    If \c NULL, the number of clusters will not be returned.
 * \param quality The quality of the partition, in terms of the objective
 *    function as included in the documentation. If \c NULL the quality will
 *    not be calculated.
 * \return Error code.
 *
 * Time complexity: near linear on sparse graphs.
 *
 * \sa \ref igraph_community_leiden_simple() for a simplified interface
 * that allows specifying an objective function directly and does not require
 * vertex weights.
 *
 * \example examples/simple/igraph_community_leiden.c
 */
igraph_error_t igraph_community_leiden(
        const igraph_t *graph,
        const igraph_vector_t *edge_weights,
        const igraph_vector_t *vertex_out_weights,
        const igraph_vector_t *vertex_in_weights,
        igraph_real_t resolution,
        igraph_real_t beta,
        igraph_bool_t start,
        igraph_int_t n_iterations,
        igraph_vector_int_t *membership,
        igraph_int_t *nb_clusters,
        igraph_real_t *quality) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    const igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_t *i_edge_weights, *i_vertex_out_weights, *i_vertex_in_weights;
    igraph_int_t i_nb_clusters;

    if (!nb_clusters) {
        nb_clusters = &i_nb_clusters;
    }

    if (start) {
        if (!membership) {
            IGRAPH_ERROR("Cannot start optimization if membership is missing.", IGRAPH_EINVAL);
        }

        if (igraph_vector_int_size(membership) != vcount) {
            IGRAPH_ERROR("Membership vector length does not equal the number of vertices.", IGRAPH_EINVAL);
        }
    } else {
        if (!membership)
            IGRAPH_ERROR("Membership vector should be supplied and initialized, "
                         "even when not starting optimization from it.", IGRAPH_EINVAL);

        IGRAPH_CHECK(igraph_vector_int_range(membership, 0, vcount));
    }

    /* Check edge weights to possibly use default. */
    if (!edge_weights) {
        i_edge_weights = IGRAPH_CALLOC(1, igraph_vector_t);
        IGRAPH_CHECK_OOM(i_edge_weights, "Leiden algorithm failed, could not allocate memory for edge weights.");
        IGRAPH_FINALLY(igraph_free, i_edge_weights);
        IGRAPH_CHECK(igraph_vector_init(i_edge_weights, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_vector_destroy, i_edge_weights);
        igraph_vector_fill(i_edge_weights, 1);
    } else {
        if (igraph_vector_size(edge_weights) != ecount) {
            IGRAPH_ERRORF("Edge weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                          IGRAPH_EINVAL, igraph_vector_size(edge_weights), ecount);
        }
        i_edge_weights = (igraph_vector_t*)edge_weights;
    }

    /* Check vertex out-weights to possibly use default. */
    if (!vertex_out_weights) {
        i_vertex_out_weights = IGRAPH_CALLOC(1, igraph_vector_t);
        IGRAPH_CHECK_OOM(i_vertex_out_weights, "Leiden algorithm failed, could not allocate memory for vertex weights.");
        IGRAPH_FINALLY(igraph_free, i_vertex_out_weights);
        IGRAPH_VECTOR_INIT_FINALLY(i_vertex_out_weights, vcount);
        igraph_vector_fill(i_vertex_out_weights, 1);
    } else {
        if (igraph_vector_size(vertex_out_weights) != vcount) {
            IGRAPH_ERRORF("Vertex %sweight vector length (%" IGRAPH_PRId ") does not match number of vertices (%" IGRAPH_PRId ").",
                          IGRAPH_EINVAL,
                          directed ? "out-" : "",
                          igraph_vector_size(vertex_out_weights), vcount);
        }
        i_vertex_out_weights = (igraph_vector_t*)vertex_out_weights;
    }

    if (directed) {
        /* When in-weights are not given for a directed graph,
         * assume that they are the same as the out-weights.
         * This effectively ignores edge directions. */
        if (vertex_in_weights) {
            if (igraph_vector_size(vertex_in_weights) != vcount) {
                IGRAPH_ERRORF("Vertex in-weight vector length (%" IGRAPH_PRId ") does not match number of vertices (%" IGRAPH_PRId ").",
                              IGRAPH_EINVAL,
                              igraph_vector_size(vertex_in_weights), vcount);
            }
            i_vertex_in_weights = (igraph_vector_t*)vertex_in_weights;
        } else {
            i_vertex_in_weights = i_vertex_out_weights;
        }
    } else {
        /* In-weights must be NULL in the undirected case. */
        if (vertex_in_weights) {
            IGRAPH_ERROR("Vertex in-weights must not be given for undirected graphs.", IGRAPH_EINVAL);
        } else {
            i_vertex_in_weights = NULL;
        }
    }

    /* Perform actual Leiden algorithm iteratively. We either
     * perform a fixed number of iterations, or we perform
     * iterations until the quality remains unchanged. Even if
     * a single iteration did not change anything, a subsequent
     * iteration may still find some improvement. This is because
     * each iteration explores different subsets of vertices.
     */
    igraph_bool_t changed = true;
    for (igraph_int_t itr = 0;
         n_iterations < 0 ? changed : itr < n_iterations;
         itr++) {
        IGRAPH_CHECK(community_leiden(graph,
                                      i_edge_weights, i_vertex_out_weights, i_vertex_in_weights,
                                      resolution, beta,
                                      membership, nb_clusters, quality, &changed));
    }

    if (!edge_weights) {
        igraph_vector_destroy(i_edge_weights);
        IGRAPH_FREE(i_edge_weights);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (!vertex_out_weights) {
        igraph_vector_destroy(i_vertex_out_weights);
        IGRAPH_FREE(i_vertex_out_weights);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_community_leiden_simple
 * \brief Finding community structure using the Leiden algorithm, simple interface.
 *
 * This is a simplified interface to \ref igraph_community_leiden() for
 * convenience purposes. Instead of requiring vertex weights, it allows
 * choosing from a set of objective functions to maximize. It implements
 * these objective functions by passing suitable vertex weights to
 * \ref igraph_community_leiden(), as explained in the documentation of
 * that function.
 *
 * \param graph The input graph. May be directed or undirected.
 * \param weights The edge weights. If \c NULL, all weights are assumed to be 1.
 * \param objective The objective function to maximize.
 *    \clist
 *    \cli IGRAPH_LEIDEN_OBJECTIVE_MODULARITY
 *      Use the generalized modularity, defined as
 *      <code>Q = 1/(2m) sum_ij (A_ij - γ k_i k_j / (2m)) δ(c_i, c_j)</code>
 *      for undirected graphs and as
 *      <code>Q = 1/m sum_ij (A_ij - γ k^out_i k^in_j / m) δ(c_i, c_j)</code>
 *      for directed graphs. This effectively uses a multigraph configuration
 *      model as the null model. Edge weights must not be negative.
 *    \cli IGRAPH_LEIDEN_OBJECTIVE_CPM
 *      Use the constant Potts model, whose objective function is defined as
 *      <code>Q = 1/(2m) sum_ij (A_ij - γ) δ(c_i, c_j)</code>
 *      for undirected graphs and as
 *      <code>Q = 1/m sum_ij (A_ij - γ) δ(c_i, c_j)</code>
 *      for directed graphs. Edge weights are allowed to be negative.
 *      Edge directions have no impact on the result.
 *    \cli IGRAPH_LEIDEN_OBJECTIVE_ER
 *      Use an objective function based on the multigraph Erdős-Rényi G(n,p)
 *      null model, defined as
 *      <code>Q = 1/(2m) sum_ij (A_ij - γ p) δ(c_i, c_j)</code>
 *      for undirected graphs and as
 *      <code>Q = 1/m sum_ij (A_ij - γ p) δ(c_i, c_j)</code>
 *      for directed graphs. \c p is the weighted density, i.e. the average
 *      link strength between all vertex pairs (whether adjacent or not).
 *      Edge weights must not be negative. Edge directions have no impact on
 *      the result.
 *    \endclist
 *    In the above formulas, \c A is the adjacency matrix, \c m is the total
 *    edge weight, \c k are the (out- and in-) degrees, \c γ is the resolution
 *    parameter, and <code>δ(c_i, c_j)</code> is 1 if vertices \c i and \c j
 *    are in the same community and 0 otherwise. Edge directions are only
 *    relevant with \c IGRAPH_LEIDEN_OBJECTIVE_MODULARITY. The other two
 *    objective functions are equivalent between directed and undirected graphs:
 *    the formal difference is due to each edge being included twice in
 *    undirected (symmetric) adjacency matrices.
 * \param resolution The resolution parameter, which is represented by γ in
 *    the objective functions detailed above.
 * \param beta The randomness used in the refinement step when merging. A small
 *    amount of randomness (\c beta = 0.01) typically works well.
 * \param start Start from membership vector. If this is true, the optimization
 *    will start from the provided membership vector. If this is false, the
 *    optimization will start from a singleton partition.
 * \param n_iterations Iterate the core Leiden algorithm the indicated number
 *    of times. If this is a negative number, it will continue iterating until
 *    an iteration did not change the clustering. Two iterations are often
 *    sufficient, thus 2 is a reasonable default.
 * \param membership The membership vector. If \p start is set to \c false,
 *    it will be resized appropriately. If \p start is \c true, it must be
 *    a valid membership vector for the given \p graph.
 * \param nb_clusters The number of clusters contained in the final \p membership.
 *    If \c NULL, the number of clusters will not be returned.
 * \param quality The quality of the partition, in terms of the objective
 *    function selected by \p objective. If \c NULL the quality will
 *    not be calculated.
 * \return Error code.
 *
 * Time complexity: near linear on sparse graphs.
 *
 * \sa \ref igraph_community_leiden() for a more flexible interface that
 * allows specifying raw vertex weights.
 */
igraph_error_t igraph_community_leiden_simple(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_leiden_objective_t objective,
        igraph_real_t resolution,
        igraph_real_t beta,
        igraph_bool_t start,
        igraph_int_t n_iterations,
        igraph_vector_int_t *membership,
        igraph_int_t *nb_clusters,
        igraph_real_t *quality) {

    const igraph_int_t vcount = igraph_vcount(graph);
    const igraph_int_t ecount = igraph_ecount(graph);
    const igraph_bool_t directed = igraph_is_directed(graph);
    igraph_vector_t vertex_out_weights, vertex_in_weights;
    igraph_vector_int_t i_membership, *p_membership;
    igraph_real_t min_weight = IGRAPH_INFINITY;

    /* Basic weight vector validation, calculate properties used for validation steps
     * specific to different objective functions. */
    if (weights) {
        if (igraph_vector_size(weights) != ecount) {
            IGRAPH_ERROR("Edge weight vector length does not match number of edges.", IGRAPH_EINVAL);
        }
        for (igraph_int_t i=0; i < ecount; i++) {
            igraph_real_t w = VECTOR(*weights)[i];
            if (w < min_weight) {
                min_weight = w;
            }
            if (! isfinite(w)) {
                IGRAPH_ERRORF("Edge weights must not be infinite or NaN, got %g.",
                              IGRAPH_EINVAL, w);
            }
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&vertex_out_weights, vcount);
    if (directed) {
        IGRAPH_VECTOR_INIT_FINALLY(&vertex_in_weights, vcount);
    }

    /* igraph_community_leiden() always requires an initialized membership vector
     * of the correct size to be given. We relax this requirement to the case
     * when start = true. */
    if (start) {
        if (!membership) {
            IGRAPH_ERROR("Requesting to start the computation from a specific "
                         "community assignment, but no membership vector given.",
                         IGRAPH_EINVAL);
        }
        if (igraph_vector_int_size(membership) != vcount) {
            IGRAPH_ERRORF("Requesting to start the computation from a specific "
                          "community assignment, but the given membership vector "
                          "has a different size (%" IGRAPH_PRId " than the vertex "
                          "count (%" IGRAPH_PRId ").",
                          IGRAPH_EINVAL,
                          igraph_vector_int_size(membership), vcount);
        }
        p_membership = membership;
    } else {
        if (!membership) {
            IGRAPH_VECTOR_INT_INIT_FINALLY(&i_membership, vcount);
            p_membership = &i_membership;
        } else {
            IGRAPH_CHECK(igraph_vector_int_resize(membership, vcount));
            p_membership = membership;
        }
    }

    switch (objective) {
    case IGRAPH_LEIDEN_OBJECTIVE_MODULARITY:
        if (min_weight < 0) {
            IGRAPH_ERRORF("Edge weights must not be negative for Leiden community "
                          "detection with modularity objective function, got %g.",
                          IGRAPH_EINVAL,
                          min_weight);
        }

        IGRAPH_CHECK(igraph_strength(
            graph, &vertex_out_weights,
            igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS, weights));
        if (directed) {
            IGRAPH_CHECK(igraph_strength(
                graph, &vertex_in_weights,
                igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS, weights));
        }

        /* If directed, the sum of vertex_out_weights is the total edge weight.
         * If undirected, it is twice the total edge weight. */
        resolution /= igraph_vector_sum(&vertex_out_weights);

        break;

    case IGRAPH_LEIDEN_OBJECTIVE_CPM:
        /* TODO: Potential minor optimization is to use the same vector for both. */
        igraph_vector_fill(&vertex_out_weights, 1);
        if (directed) {
            igraph_vector_fill(&vertex_in_weights, 1);
        }

        break;

    case IGRAPH_LEIDEN_OBJECTIVE_ER:
        if (min_weight < 0) {
            IGRAPH_ERRORF("Edge weights must not be negative for Leiden community "
                          "detection with ER objective function, got %g.",
                          IGRAPH_EINVAL,
                          min_weight);
        }

        /* TODO: Potential minor optimization is to use the same vector for both. */
        igraph_vector_fill(&vertex_out_weights, 1);
        if (directed) {
            igraph_vector_fill(&vertex_in_weights, 1);
        }

        {
            igraph_real_t p;
            /* Note: Loops must be allowed, as the aggregation step of the
             * algorithm effectively creates them. */
            IGRAPH_CHECK(igraph_density(graph, weights, &p, /* loops */ true));
            resolution *= p;
        }

        break;


    default:
        IGRAPH_ERROR("Invalid objective function for Leiden community detection.",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_community_leiden(
        graph, weights,
        &vertex_out_weights, directed ? &vertex_in_weights : NULL,
        resolution, beta, start, n_iterations, p_membership, nb_clusters, quality));

    if (!membership) {
        igraph_vector_int_destroy(&i_membership);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (directed) {
        igraph_vector_destroy(&vertex_in_weights);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&vertex_out_weights);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
