/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_stack.h"
#include "igraph_vector.h"
#include "igraph_constructors.h"

#include "core/interruption.h"

/* Move nodes in order to improve the quality of a partition.
 *
 * This function considers each node and greedily moves it to a neighboring
 * community that maximizes the improvement in the quality of a partition.
 *
 * The nodes are examined in a queue, and initially all nodes are put in the
 * queue in a random order. Nodes are popped from the queue when they are
 * examined, and only neighbors of nodes that are moved (which are not part of
 * the cluster the node was moved to) are pushed to the queue again.
 *
 * The \c membership vector is used as the starting point to move around nodes,
 * and is updated in-place.
 *
 */
static int igraph_i_community_leiden_fastmovenodes(
        const igraph_t *graph,
        const igraph_inclist_t *edges_per_node,
        const igraph_vector_t *edge_weights, const igraph_vector_t *node_weights,
        const igraph_real_t resolution_parameter,
        igraph_integer_t *nb_clusters,
        igraph_vector_t *membership) {

    igraph_dqueue_t unstable_nodes;
    igraph_real_t max_diff = 0.0, diff = 0.0;
    igraph_integer_t n = igraph_vcount(graph);
    igraph_vector_bool_t neighbor_cluster_added, node_is_stable;
    igraph_vector_t node_order, cluster_weights, edge_weights_per_cluster, neighbor_clusters;
    igraph_vector_int_t nb_nodes_per_cluster;
    igraph_stack_t empty_clusters;
    long int i, j, c, nb_neigh_clusters;

    /* Initialize queue of unstable nodes and whether node is stable. Only
     * unstable nodes are in the queue. */
    IGRAPH_CHECK(igraph_vector_bool_init(&node_is_stable, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &node_is_stable);

    IGRAPH_CHECK(igraph_dqueue_init(&unstable_nodes, n));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &unstable_nodes);

    /* Shuffle nodes */
    IGRAPH_CHECK(igraph_vector_init_seq(&node_order, 0, n - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &node_order);
    IGRAPH_CHECK(igraph_vector_shuffle(&node_order));

    /* Add to the queue */
    for (i = 0; i < n; i++) {
        igraph_dqueue_push(&unstable_nodes, (long int)VECTOR(node_order)[i]);
    }

    /* Initialize cluster weights and nb nodes */
    IGRAPH_CHECK(igraph_vector_init(&cluster_weights, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &cluster_weights);
    IGRAPH_CHECK(igraph_vector_int_init(&nb_nodes_per_cluster, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &nb_nodes_per_cluster);
    for (i = 0; i < n; i++) {
        c = (long int)VECTOR(*membership)[i];
        VECTOR(cluster_weights)[c] += VECTOR(*node_weights)[i];
        VECTOR(nb_nodes_per_cluster)[c] += 1;
    }

    /* Initialize empty clusters */
    IGRAPH_CHECK(igraph_stack_init(&empty_clusters, n));
    IGRAPH_FINALLY(igraph_stack_destroy, &empty_clusters);
    for (c = 0; c < n; c++)
        if (VECTOR(nb_nodes_per_cluster)[c] == 0) {
            igraph_stack_push(&empty_clusters, c);
        }

    /* Initialize vectors to be used in calculating differences */
    IGRAPH_CHECK(igraph_vector_init(&edge_weights_per_cluster, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &edge_weights_per_cluster);

    /* Initialize neighboring cluster */
    IGRAPH_CHECK(igraph_vector_bool_init(&neighbor_cluster_added, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &neighbor_cluster_added);
    IGRAPH_CHECK(igraph_vector_init(&neighbor_clusters, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &neighbor_clusters);

    /* Iterate while the queue is not empty */
    j = 0;
    while (!igraph_dqueue_empty(&unstable_nodes)) {
        long int v = (long int) igraph_dqueue_pop(&unstable_nodes);
        long int best_cluster, current_cluster = VECTOR(*membership)[v];
        long int degree, i;
        igraph_vector_int_t *edges;

        /* Remove node from current cluster */
        VECTOR(cluster_weights)[current_cluster] -= VECTOR(*node_weights)[v];
        VECTOR(nb_nodes_per_cluster)[current_cluster]--;
        if (VECTOR(nb_nodes_per_cluster)[current_cluster] == 0) {
            IGRAPH_CHECK(igraph_stack_push(&empty_clusters, current_cluster));
        }

        /* Find out neighboring clusters */
        c = (long int) igraph_stack_top(&empty_clusters);
        VECTOR(neighbor_clusters)[0] = c;
        VECTOR(neighbor_cluster_added)[c] = 1;
        nb_neigh_clusters = 1;

        /* Determine the edge weight to each neighboring cluster */
        edges = igraph_inclist_get(edges_per_node, v);
        degree = igraph_vector_int_size(edges);
        for (i = 0; i < degree; i++) {
            long int e = VECTOR(*edges)[i];
            long int u = (long int)IGRAPH_OTHER(graph, e, v);
            if (u != v) {
                c = VECTOR(*membership)[u];
                if (!VECTOR(neighbor_cluster_added)[c]) {
                    VECTOR(neighbor_cluster_added)[c] = 1;
                    VECTOR(neighbor_clusters)[nb_neigh_clusters++] = c;
                }
                VECTOR(edge_weights_per_cluster)[c] += VECTOR(*edge_weights)[e];
            }
        }

        /* Calculate maximum diff */
        best_cluster = current_cluster;
        max_diff = VECTOR(edge_weights_per_cluster)[current_cluster] - VECTOR(*node_weights)[v] * VECTOR(cluster_weights)[current_cluster] * resolution_parameter;
        for (i = 0; i < nb_neigh_clusters; i++) {
            c = VECTOR(neighbor_clusters)[i];
            diff = VECTOR(edge_weights_per_cluster)[c] - VECTOR(*node_weights)[v] * VECTOR(cluster_weights)[c] * resolution_parameter;
            if (diff > max_diff) {
                best_cluster = c;
                max_diff = diff;
            }
            VECTOR(edge_weights_per_cluster)[c] = 0.0;
            VECTOR(neighbor_cluster_added)[c] = 0;
        }

        /* Move node to best cluster */
        VECTOR(cluster_weights)[best_cluster] += VECTOR(*node_weights)[v];
        VECTOR(nb_nodes_per_cluster)[best_cluster]++;
        if (best_cluster == igraph_stack_top(&empty_clusters)) {
            igraph_stack_pop(&empty_clusters);
        }

        /* Mark node as stable */
        VECTOR(node_is_stable)[v] = 1;

        /* Add stable neighbours that are not part of the new cluster to the queue */
        if (best_cluster != current_cluster) {
            VECTOR(*membership)[v] = best_cluster;

            for (i = 0; i < degree; i++) {
                long int e = VECTOR(*edges)[i];
                long int u = (long int) IGRAPH_OTHER(graph, e, v);
                if (VECTOR(node_is_stable)[u] && VECTOR(*membership)[u] != best_cluster) {
                    IGRAPH_CHECK(igraph_dqueue_push(&unstable_nodes, u));
                    VECTOR(node_is_stable)[u] = 0;
                }
            }
        }

        j++;
        if (j > 10000) {
            IGRAPH_ALLOW_INTERRUPTION();
            j = 0;
        }
    }

    IGRAPH_CHECK(igraph_reindex_membership(membership, NULL, nb_clusters));

    igraph_vector_destroy(&neighbor_clusters);
    igraph_vector_bool_destroy(&neighbor_cluster_added);
    igraph_vector_destroy(&edge_weights_per_cluster);
    igraph_stack_destroy(&empty_clusters);
    igraph_vector_int_destroy(&nb_nodes_per_cluster);
    igraph_vector_destroy(&cluster_weights);
    igraph_vector_destroy(&node_order);
    igraph_dqueue_destroy(&unstable_nodes);
    igraph_vector_bool_destroy(&node_is_stable);

    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

/* Clean a refined membership vector.
 *
 * This function examines all nodes in \c node_subset and updates \c
 * refined_membership to ensure that the clusters are numbered consecutively,
 * starting from \c nb_refined_clusters. The \c nb_refined_clusters is also
 * updated itself. If C is the initial \c nb_refined_clusters and C' the
 * resulting \c nb_refined_clusters, then nodes in \c node_subset are numbered
 * C, C + 1, ..., C' - 1.
 */
static int igraph_i_community_leiden_clean_refined_membership(
        const igraph_vector_t* node_subset,
        igraph_vector_t *refined_membership,
        igraph_integer_t* nb_refined_clusters) {
    long int i, n = igraph_vector_size(node_subset);
    igraph_vector_t new_cluster;

    IGRAPH_CHECK(igraph_vector_init(&new_cluster, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &new_cluster);

    /* Clean clusters. We will store the new cluster + 1 so that cluster == 0
     * indicates that no membership was assigned yet. */
    *nb_refined_clusters += 1;
    for (i = 0; i < n; i++) {
        long int v = (long int) VECTOR(*node_subset)[i];
        long int c = (long int) VECTOR(*refined_membership)[v];
        if (VECTOR(new_cluster)[c] == 0) {
            VECTOR(new_cluster)[c] = (igraph_real_t)(*nb_refined_clusters);
            *nb_refined_clusters += 1;
        }
    }

    /* Assign new cluster */
    for (i = 0; i < n; i++) {
        long int v = (long int) VECTOR(*node_subset)[i];
        long int c = (long int) VECTOR(*refined_membership)[v];
        VECTOR(*refined_membership)[v] = VECTOR(new_cluster)[c] - 1;
    }
    /* We used the cluster + 1, so correct */
    *nb_refined_clusters -= 1;

    igraph_vector_destroy(&new_cluster);

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Merge nodes for a subset of the nodes. This is used to refine a partition.
 *
 * The nodes included in \c node_subset are assumed to be the nodes i for which
 * membership[i] = cluster_subset.
 *
 * All nodes in \c node_subset are initialized to a singleton partition in \c
 * refined_membership. Only singleton clusters can be merged if they are
 * sufficiently well connected to the current subgraph induced by \c
 * node_subset.
 *
 * We only examine each node once. Instead of greedily choosing the maximum
 * possible cluster to merge with, the cluster is chosen randomly among all
 * possibilities that do not decrease the quality of the partition. The
 * probability of choosing a certain cluster is proportional to exp(diff/beta).
 * For beta to 0 this converges to selecting a cluster with the maximum
 * improvement. For beta to infinity this converges to a uniform distribution
 * among all eligible clusters.
 *
 * The \c refined_membership is updated for node in \c node_subset. The number
 * of refined clusters, \c nb_refined_clusters is used to set the actual refined
 * cluster membership and is updated after this routine. Within each cluster
 * (i.e. for a given \c node_subset), the refined membership is initially simply
 * set to 0, ..., n - 1 (for n nodes in \c node_subset). However, for each \c
 * node_subset the refined membership should of course be unique. Hence, after
 * merging, the refined membership starts with \c nb_refined_clusters, which is
 * also updated to ensure that the resulting \c nb_refined_clusters counts all
 * refined clusters that have already been processed. See
 * igraph_i_community_leiden_clean_refined_membership for more information about
 * this aspect.
 */
static int igraph_i_community_leiden_mergenodes(
        const igraph_t *graph,
        const igraph_inclist_t *edges_per_node,
        const igraph_vector_t *edge_weights, const igraph_vector_t *node_weights,
        const igraph_vector_t *node_subset,
        const igraph_vector_t *membership,
        const igraph_integer_t cluster_subset,
        const igraph_real_t resolution_parameter,
        const igraph_real_t beta,
        igraph_integer_t *nb_refined_clusters,
        igraph_vector_t *refined_membership) {
    igraph_vector_t node_order;
    igraph_vector_bool_t non_singleton_cluster, neighbor_cluster_added;
    igraph_real_t max_diff, total_cum_trans_diff, diff = 0.0, total_node_weight = 0.0;
    igraph_integer_t n = igraph_vector_size(node_subset);
    igraph_vector_t cluster_weights, cum_trans_diff, edge_weights_per_cluster, external_edge_weight_per_cluster_in_subset, neighbor_clusters;
    igraph_vector_int_t *edges, nb_nodes_per_cluster;
    long int i, j, degree, nb_neigh_clusters;

    /* Initialize cluster weights */
    IGRAPH_CHECK(igraph_vector_init(&cluster_weights, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &cluster_weights);

    /* Initialize number of nodes per cluster */
    IGRAPH_CHECK(igraph_vector_int_init(&nb_nodes_per_cluster, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &nb_nodes_per_cluster);

    /* Initialize external edge weight per cluster in subset */
    IGRAPH_CHECK(igraph_vector_init(&external_edge_weight_per_cluster_in_subset, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &external_edge_weight_per_cluster_in_subset);

    /* Initialize administration for a singleton partition */
    for (i = 0; i < n; i++) {
        long int v = (long int) VECTOR(*node_subset)[i];
        VECTOR(*refined_membership)[v] = i;
        VECTOR(cluster_weights)[i] += VECTOR(*node_weights)[v];
        VECTOR(nb_nodes_per_cluster)[i] += 1;
        total_node_weight += VECTOR(*node_weights)[v];

        /* Find out neighboring clusters */
        edges = igraph_inclist_get(edges_per_node, v);
        degree = igraph_vector_int_size(edges);
        for (j = 0; j < degree; j++) {
            long int e = VECTOR(*edges)[j];
            long int u = (long int)IGRAPH_OTHER(graph, e, v);
            if (u != v && VECTOR(*membership)[u] == cluster_subset) {
                VECTOR(external_edge_weight_per_cluster_in_subset)[i] += VECTOR(*edge_weights)[e];
            }
        }
    }

    /* Shuffle nodes */
    IGRAPH_CHECK(igraph_vector_copy(&node_order, node_subset));
    IGRAPH_FINALLY(igraph_vector_destroy, &node_order);
    IGRAPH_CHECK(igraph_vector_shuffle(&node_order));

    /* Initialize non singleton clusters */
    IGRAPH_CHECK(igraph_vector_bool_init(&non_singleton_cluster, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &non_singleton_cluster);

    /* Initialize vectors to be used in calculating differences */
    IGRAPH_CHECK(igraph_vector_init(&edge_weights_per_cluster, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &edge_weights_per_cluster);

    /* Initialize neighboring cluster */
    IGRAPH_CHECK(igraph_vector_bool_init(&neighbor_cluster_added, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &neighbor_cluster_added);
    IGRAPH_CHECK(igraph_vector_init(&neighbor_clusters, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &neighbor_clusters);

    /* Initialize cumulative transformed difference */
    IGRAPH_CHECK(igraph_vector_init(&cum_trans_diff, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &cum_trans_diff);

    RNG_BEGIN();

    for (i = 0; i < n; i++) {
        long int v = (long int) VECTOR(node_order)[i];
        long int chosen_cluster, best_cluster, current_cluster = (long int) VECTOR(*refined_membership)[v];

        if (!VECTOR(non_singleton_cluster)[current_cluster] &&
            (VECTOR(external_edge_weight_per_cluster_in_subset)[current_cluster] >=
             VECTOR(cluster_weights)[current_cluster] * (total_node_weight - VECTOR(cluster_weights)[current_cluster]) * resolution_parameter)) {
            /* Remove node from current cluster, which is then a singleton by
             * definition. */
            VECTOR(cluster_weights)[current_cluster] = 0.0;
            VECTOR(nb_nodes_per_cluster)[current_cluster] = 0;

            /* Find out neighboring clusters */
            edges = igraph_inclist_get(edges_per_node, v);
            degree = igraph_vector_int_size(edges);

            /* Also add current cluster to ensure it can be chosen. */
            VECTOR(neighbor_clusters)[0] = current_cluster;
            VECTOR(neighbor_cluster_added)[current_cluster] = 1;
            nb_neigh_clusters = 1;
            for (j = 0; j < degree; j++) {
                long int e = (long int)VECTOR(*edges)[j];
                long int u = (long int)IGRAPH_OTHER(graph, e, v);
                if (u != v && VECTOR(*membership)[u] == cluster_subset) {
                    long int c = VECTOR(*refined_membership)[u];
                    if (!VECTOR(neighbor_cluster_added)[c]) {
                        VECTOR(neighbor_cluster_added)[c] = 1;
                        VECTOR(neighbor_clusters)[nb_neigh_clusters++] = c;
                    }
                    VECTOR(edge_weights_per_cluster)[c] += VECTOR(*edge_weights)[e];
                }
            }

            /* Calculate diffs */
            best_cluster = current_cluster;
            max_diff = 0.0;
            total_cum_trans_diff = 0.0;
            for (j = 0; j < nb_neigh_clusters; j++) {
                long int c = (long int) VECTOR(neighbor_clusters)[j];
                if (VECTOR(external_edge_weight_per_cluster_in_subset)[c] >= VECTOR(cluster_weights)[c] * (total_node_weight - VECTOR(cluster_weights)[c]) * resolution_parameter) {
                    diff = VECTOR(edge_weights_per_cluster)[c] - VECTOR(*node_weights)[v] * VECTOR(cluster_weights)[c] * resolution_parameter;

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
                VECTOR(neighbor_cluster_added)[c] = 0;
            }

            /* Determine the neighboring cluster to which the currently selected node
             * will be moved.
             */
            if (total_cum_trans_diff < IGRAPH_INFINITY) {
                igraph_real_t r = RNG_UNIF(0, total_cum_trans_diff);
                long int chosen_idx;
                igraph_vector_binsearch_slice(&cum_trans_diff, r, &chosen_idx, 0, nb_neigh_clusters);
                chosen_cluster = VECTOR(neighbor_clusters)[chosen_idx];
            } else {
                chosen_cluster = best_cluster;
            }

            /* Move node to randomly chosen cluster */
            VECTOR(cluster_weights)[chosen_cluster] += VECTOR(*node_weights)[v];
            VECTOR(nb_nodes_per_cluster)[chosen_cluster]++;

            for (j = 0; j < degree; j++) {
                long int e = (long int) VECTOR(*edges)[j];
                long int u = (long int) IGRAPH_OTHER(graph, e, v);
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

                VECTOR(non_singleton_cluster)[chosen_cluster] = 1;
            }
        } /* end if singleton and may be merged */
    }

    RNG_END();

    IGRAPH_CHECK(igraph_i_community_leiden_clean_refined_membership(node_subset, refined_membership, nb_refined_clusters));

    igraph_vector_destroy(&cum_trans_diff);
    igraph_vector_destroy(&neighbor_clusters);
    igraph_vector_bool_destroy(&neighbor_cluster_added);
    igraph_vector_destroy(&edge_weights_per_cluster);
    igraph_vector_bool_destroy(&non_singleton_cluster);
    igraph_vector_destroy(&node_order);
    igraph_vector_destroy(&external_edge_weight_per_cluster_in_subset);
    igraph_vector_int_destroy(&nb_nodes_per_cluster);
    igraph_vector_destroy(&cluster_weights);

    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

/* Create clusters out of a membership vector.
 *
 * The cluster pointer vector should be initialized for all entries of the
 * membership vector, no range checking is performed. If a vector for a cluster
 * does not yet exist it will be created and initialized. If a vector for a
 * cluster already does exist it will not be emptied on first use. Hence, it
 * should be ensured that all clusters are always properly empty (or
 * non-existing) before calling this function.
 */
static int igraph_i_community_get_clusters(const igraph_vector_t *membership, igraph_vector_ptr_t *clusters) {
    long int i, c, n = igraph_vector_size(membership);
    igraph_vector_t *cluster;
    for (i = 0; i < n; i++) {
        /* Get cluster for node i */
        c = VECTOR(*membership)[i];
        cluster = (igraph_vector_t*)VECTOR(*clusters)[c];

        /* No cluster vector exists yet, so we create a new one */
        if (!cluster) {
            cluster = IGRAPH_CALLOC(1, igraph_vector_t);
            if (cluster == 0) {
                IGRAPH_ERROR("Cannot allocate memory for assigning cluster", IGRAPH_ENOMEM);
            }
            IGRAPH_CHECK(igraph_vector_init(cluster, 0));
            VECTOR(*clusters)[c] = cluster;
        }

        /* Add node i to cluster vector */
        IGRAPH_CHECK(igraph_vector_push_back(cluster, i));
    }

    return IGRAPH_SUCCESS;
}

/* Aggregate the graph based on the \c refined membership while setting the
 * membership of each aggregated node according to the \c membership.
 *
 * Technically speaking we have that
 * aggregated_membership[refined_membership[v]] = membership[v] for each node v.
 *
 * The new aggregated graph is returned in \c aggregated_graph. This graph
 * object should not yet be initialized, `igraph_create` is called on it, and
 * responsibility for destroying the object lies with the calling method
 *
 * The remaining results, aggregated_edge_weights, aggregate_node_weights and
 * aggregated_membership are all expected to be initialized.
 *
 */
static int igraph_i_community_leiden_aggregate(
    const igraph_t *graph, const igraph_inclist_t *edges_per_node, const igraph_vector_t *edge_weights, const igraph_vector_t *node_weights,
    const igraph_vector_t *membership, const igraph_vector_t *refined_membership, const igraph_integer_t nb_refined_clusters,
    igraph_t *aggregated_graph, igraph_vector_t *aggregated_edge_weights, igraph_vector_t *aggregated_node_weights, igraph_vector_t *aggregated_membership) {
    igraph_vector_t aggregated_edges, edge_weight_to_cluster;
    igraph_vector_ptr_t refined_clusters;
    igraph_vector_int_t *incident_edges;
    igraph_vector_t neighbor_clusters;
    igraph_vector_bool_t neighbor_cluster_added;
    long int i, j, c, degree, nb_neigh_clusters;

    /* Get refined clusters */
    IGRAPH_CHECK(igraph_vector_ptr_init(&refined_clusters, nb_refined_clusters));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&refined_clusters, igraph_vector_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &refined_clusters);
    IGRAPH_CHECK(igraph_i_community_get_clusters(refined_membership, &refined_clusters));

    /* Initialize new edges */
    IGRAPH_CHECK(igraph_vector_init(&aggregated_edges, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &aggregated_edges);

    /* We clear the aggregated edge weights, we will push each new edge weight */
    igraph_vector_clear(aggregated_edge_weights);
    /* Simply resize the aggregated node weights and membership, they can be set
     * directly */
    IGRAPH_CHECK(igraph_vector_resize(aggregated_node_weights, nb_refined_clusters));
    IGRAPH_CHECK(igraph_vector_resize(aggregated_membership, nb_refined_clusters));

    IGRAPH_CHECK(igraph_vector_init(&edge_weight_to_cluster, nb_refined_clusters));
    IGRAPH_FINALLY(igraph_vector_destroy, &edge_weight_to_cluster);

    /* Initialize neighboring cluster */
    IGRAPH_CHECK(igraph_vector_bool_init(&neighbor_cluster_added, nb_refined_clusters));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &neighbor_cluster_added);
    IGRAPH_CHECK(igraph_vector_init(&neighbor_clusters, nb_refined_clusters));
    IGRAPH_FINALLY(igraph_vector_destroy, &neighbor_clusters);

    /* Check per cluster */
    for (c = 0; c < nb_refined_clusters; c++) {
        igraph_vector_t* refined_cluster = (igraph_vector_t*)VECTOR(refined_clusters)[c];
        long int n_c = igraph_vector_size(refined_cluster);
        long int v = -1;

        /* Calculate the total edge weight to other clusters */
        VECTOR(*aggregated_node_weights)[c] = 0.0;
        nb_neigh_clusters = 0;
        for (i = 0; i < n_c; i++) {
            v = (long int) VECTOR(*refined_cluster)[i];
            incident_edges = igraph_inclist_get(edges_per_node, v);
            degree = igraph_vector_int_size(incident_edges);

            for (j = 0; j < degree; j++) {
                long int e = VECTOR(*incident_edges)[j];
                long int u = (long int) IGRAPH_OTHER(graph, e, v);
                long int c2 = VECTOR(*refined_membership)[u];

                if (c2 > c) {
                    if (!VECTOR(neighbor_cluster_added)[c2]) {
                        VECTOR(neighbor_cluster_added)[c2] = 1;
                        VECTOR(neighbor_clusters)[nb_neigh_clusters++] = c2;
                    }
                    VECTOR(edge_weight_to_cluster)[c2] += VECTOR(*edge_weights)[e];
                }
            }

            VECTOR(*aggregated_node_weights)[c] += VECTOR(*node_weights)[v];
        }

        /* Add actual edges from this cluster to the other clusters */
        for (i = 0; i < nb_neigh_clusters; i++) {
            long int c2 = VECTOR(neighbor_clusters)[i];

            /* Add edge */
            IGRAPH_CHECK(igraph_vector_push_back(&aggregated_edges, c));
            IGRAPH_CHECK(igraph_vector_push_back(&aggregated_edges, c2));

            /* Add edge weight */
            IGRAPH_CHECK(igraph_vector_push_back(aggregated_edge_weights, VECTOR(edge_weight_to_cluster)[c2]));

            VECTOR(edge_weight_to_cluster)[c2] = 0.0;
            VECTOR(neighbor_cluster_added)[c2] = 0;
        }

        VECTOR(*aggregated_membership)[c] = VECTOR(*membership)[v];

    }

    igraph_vector_destroy(&neighbor_clusters);
    igraph_vector_bool_destroy(&neighbor_cluster_added);
    igraph_vector_destroy(&edge_weight_to_cluster);
    igraph_vector_ptr_destroy_all(&refined_clusters);

    IGRAPH_FINALLY_CLEAN(4);

    igraph_destroy(aggregated_graph);
    IGRAPH_CHECK(igraph_create(aggregated_graph, &aggregated_edges, nb_refined_clusters,
                               IGRAPH_UNDIRECTED));

    igraph_vector_destroy(&aggregated_edges);

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* Calculate the quality of the partition.
 *
 * The quality is defined as
 *
 * 1 / 2m sum_ij (A_ij - gamma n_i n_j)d(s_i, s_j)
 *
 * where m is the total edge weight, A_ij is the weight of edge (i, j), gamma is
 * the so-called resolution parameter, n_i is the node weight of node i, s_i is
 * the cluster of node i and d(x, y) = 1 if and only if x = y and 0 otherwise.
 *
 * Note that by setting n_i = k_i the degree of node i and dividing gamma by 2m,
 * we effectively optimize modularity. By setting n_i = 1 we optimize the
 * Constant Potts Model.
 *
 * This can be represented as a sum over clusters as
 *
 * 1 / 2m sum_c (e_c - gamma N_c^2)
 *
 * where e_c = sum_ij A_ij d(s_i, c)d(s_j, c) is (twice) the internal edge
 * weight in cluster c and N_c = sum_i n_i d(s_i, c) is the sum of the node
 * weights inside cluster c. This is how the quality is calculated in practice.
 *
 */
static int igraph_i_community_leiden_quality(
        const igraph_t *graph, const igraph_vector_t *edge_weights, const igraph_vector_t *node_weights,
        const igraph_vector_t *membership, const igraph_integer_t nb_comms, const igraph_real_t resolution_parameter,
        igraph_real_t *quality) {
    igraph_vector_t cluster_weights;
    igraph_real_t total_edge_weight = 0.0;
    igraph_eit_t eit;
    long int i, c, n = igraph_vcount(graph);;

    *quality = 0.0;

    /* Create the edgelist */
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t e = IGRAPH_EIT_GET(eit), from, to;
        IGRAPH_CHECK(igraph_edge(graph, e, &from, &to));
        total_edge_weight += VECTOR(*edge_weights)[e];
        /* We add the internal edge weights */
        if (VECTOR(*membership)[(long int) from] == VECTOR(*membership)[(long int) to]) {
            *quality += 2 * VECTOR(*edge_weights)[e];
        }
        IGRAPH_EIT_NEXT(eit);
    }
    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    /* Initialize cluster weights and nb nodes */
    IGRAPH_CHECK(igraph_vector_init(&cluster_weights, n));
    IGRAPH_FINALLY(igraph_vector_destroy, &cluster_weights);
    for (i = 0; i < n; i++) {
        c = VECTOR(*membership)[i];
        VECTOR(cluster_weights)[c] += VECTOR(*node_weights)[i];
    }

    /* We subtract gamma * N_c^2 */
    for (c = 0; c < nb_comms; c++) {
        *quality -= resolution_parameter * VECTOR(cluster_weights)[c] * VECTOR(cluster_weights)[c];
    }

    igraph_vector_destroy(&cluster_weights);
    IGRAPH_FINALLY_CLEAN(1);

    /* We normalise by 2m */
    *quality /= (2.0 * total_edge_weight);

    return IGRAPH_SUCCESS;
}

/* This is the core of the Leiden algorithm and relies on subroutines to
 * perform the three different phases: (1) local moving of nodes, (2)
 * refinement of the partition and (3) aggregation of the network based on the
 * refined partition, using the non-refined partition to create an initial
 * partition for the aggregate network.
 */
static int igraph_i_community_leiden(
        const igraph_t *graph,
        igraph_vector_t *edge_weights, igraph_vector_t *node_weights,
        const igraph_real_t resolution_parameter, const igraph_real_t beta,
        igraph_vector_t *membership, igraph_integer_t *nb_clusters, igraph_real_t *quality) {
    igraph_integer_t nb_refined_clusters;
    long int i, c, n = igraph_vcount(graph);
    igraph_t aggregated_graph, *i_graph;
    igraph_vector_t aggregated_edge_weights, aggregated_node_weights, aggregated_membership;
    igraph_vector_t *i_edge_weights, *i_node_weights, *i_membership;
    igraph_vector_t tmp_edge_weights, tmp_node_weights, tmp_membership;
    igraph_vector_t refined_membership;
    igraph_vector_int_t aggregate_node;
    igraph_vector_ptr_t clusters;
    igraph_inclist_t edges_per_node;
    igraph_bool_t continue_clustering;
    igraph_integer_t level = 0;

    /* Initialize temporary weights and membership to be used in aggregation */
    IGRAPH_CHECK(igraph_vector_init(&tmp_edge_weights, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp_edge_weights);
    IGRAPH_CHECK(igraph_vector_init(&tmp_node_weights, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp_node_weights);
    IGRAPH_CHECK(igraph_vector_init(&tmp_membership, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &tmp_membership);

    /* Initialize clusters */
    IGRAPH_CHECK(igraph_vector_ptr_init(&clusters, n));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&clusters, igraph_vector_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &clusters);

    /* Initialize aggregate nodes, which initially is identical to simply the
     * nodes in the graph. */
    IGRAPH_CHECK(igraph_vector_int_init(&aggregate_node, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &aggregate_node);
    for (i = 0; i < n; i++) {
        VECTOR(aggregate_node)[i] = i;
    }

    /* Initialize refined membership */
    IGRAPH_CHECK(igraph_vector_init(&refined_membership, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &refined_membership);

    /* Initialize aggregated graph */
    IGRAPH_CHECK(igraph_empty(&aggregated_graph, 0, IGRAPH_UNDIRECTED));
    IGRAPH_FINALLY(igraph_destroy, &aggregated_graph);

    /* Initialize aggregated edge weights */
    IGRAPH_CHECK(igraph_vector_init(&aggregated_edge_weights, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &aggregated_edge_weights);

    /* Initialize aggregated node weights */
    IGRAPH_CHECK(igraph_vector_init(&aggregated_node_weights, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &aggregated_node_weights);

    /* Initialize aggregated membership */
    IGRAPH_CHECK(igraph_vector_init(&aggregated_membership, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &aggregated_membership);

    /* Set actual graph, weights and membership to be used. */
    i_graph = (igraph_t*)graph;
    i_edge_weights = edge_weights;
    i_node_weights = node_weights;
    i_membership = membership;

    /* Clean membership and count number of *clusters */

    IGRAPH_CHECK(igraph_reindex_membership(i_membership, NULL, nb_clusters));

    if (*nb_clusters > n) {
        IGRAPH_ERROR("Too many communities in membership vector", IGRAPH_EINVAL);
    }

    do {

        /* Get incidence list for fast iteration */
        IGRAPH_CHECK(igraph_inclist_init( i_graph, &edges_per_node, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &edges_per_node);

        /* Move around the nodes in order to increase the quality */
        IGRAPH_CHECK(igraph_i_community_leiden_fastmovenodes(i_graph,
                     &edges_per_node,
                     i_edge_weights, i_node_weights,
                     resolution_parameter,
                     nb_clusters,
                     i_membership));

        /* We only continue clustering if not all clusters are represented by a
         * single node yet
         */
        continue_clustering = (*nb_clusters < igraph_vcount(i_graph));

        if (continue_clustering) {
            /* Set original membership */
            if (level > 0) {
                for (i = 0; i < n; i++) {
                    long int v_aggregate = VECTOR(aggregate_node)[i];
                    VECTOR(*membership)[i] = VECTOR(*i_membership)[v_aggregate];
                }
            }

            /* Get node sets for each cluster. */
            IGRAPH_CHECK(igraph_i_community_get_clusters(i_membership, &clusters));

            /* Ensure refined membership is correct size */
            IGRAPH_CHECK(igraph_vector_resize(&refined_membership, igraph_vcount(i_graph)));

            /* Refine each cluster */
            nb_refined_clusters = 0;
            for (c = 0; c < *nb_clusters; c++) {
                igraph_vector_t* cluster = (igraph_vector_t*)VECTOR(clusters)[c];
                IGRAPH_CHECK(igraph_i_community_leiden_mergenodes(i_graph,
                             &edges_per_node,
                             i_edge_weights, i_node_weights,
                             cluster, i_membership, c,
                             resolution_parameter, beta,
                             &nb_refined_clusters, &refined_membership));
                /* Empty cluster */
                igraph_vector_clear(cluster);
            }

            /* If refinement didn't aggregate anything, we aggregate on the basis of
             * the actual clustering */
            if (nb_refined_clusters >= igraph_vcount(i_graph)) {
                igraph_vector_update(&refined_membership, i_membership);
                nb_refined_clusters = *nb_clusters;
            }

            /* Keep track of aggregate node. */
            for (i = 0; i < n; i++) {
                /* Current aggregate node */
                igraph_integer_t v_aggregate = VECTOR(aggregate_node)[i];
                /* New aggregate node */
                VECTOR(aggregate_node)[i] = (igraph_integer_t)VECTOR(refined_membership)[v_aggregate];
            }

            IGRAPH_CHECK(igraph_i_community_leiden_aggregate(
                             i_graph, &edges_per_node, i_edge_weights, i_node_weights,
                             i_membership, &refined_membership, nb_refined_clusters,
                             &aggregated_graph, &tmp_edge_weights, &tmp_node_weights, &tmp_membership));

            /* On the lowest level, the actual graph and node and edge weights and
             * membership are used. On higher levels, we will use the aggregated graph
             * and associated vectors.
             */
            if (level == 0) {
                /* Set actual graph, weights and membership to be used. */
                i_graph = &aggregated_graph;
                i_edge_weights = &aggregated_edge_weights;
                i_node_weights = &aggregated_node_weights;
                i_membership = &aggregated_membership;
            }

            /* Update the aggregated administration. */
            IGRAPH_CHECK(igraph_vector_update(i_edge_weights, &tmp_edge_weights));
            IGRAPH_CHECK(igraph_vector_update(i_node_weights, &tmp_node_weights));
            IGRAPH_CHECK(igraph_vector_update(i_membership, &tmp_membership));

            level += 1;
        }

        /* We are done iterating, so we destroy the incidence list */
        igraph_inclist_destroy(&edges_per_node);
        IGRAPH_FINALLY_CLEAN(1);
    } while (continue_clustering);

    /* Free aggregated graph and associated vectors */
    igraph_vector_destroy(&aggregated_membership);
    igraph_vector_destroy(&aggregated_node_weights);
    igraph_vector_destroy(&aggregated_edge_weights);
    igraph_destroy(&aggregated_graph);
    IGRAPH_FINALLY_CLEAN(4);

    /* Free remaining memory */
    igraph_vector_destroy(&refined_membership);
    igraph_vector_int_destroy(&aggregate_node);
    igraph_vector_ptr_destroy_all(&clusters);
    igraph_vector_destroy(&tmp_membership);
    igraph_vector_destroy(&tmp_node_weights);
    igraph_vector_destroy(&tmp_edge_weights);
    IGRAPH_FINALLY_CLEAN(6);

    /* Calculate quality */
    if (quality) {
        IGRAPH_CHECK(igraph_i_community_leiden_quality(graph, edge_weights, node_weights, membership, *nb_clusters, resolution_parameter, quality));
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup communities
 * \function igraph_community_leiden
 * \brief Finding community structure using the Leiden algorithm.
 *
 * This function implements the Leiden algorithm for finding community
 * structure, see Traag, V. A., Waltman, L., &amp; van Eck, N. J. (2019). From
 * Louvain to Leiden: guaranteeing well-connected communities. Scientific
 * reports, 9(1), 5233. http://dx.doi.org/10.1038/s41598-019-41695-z
 *
 * </para><para>
 * It is similar to the multilevel algorithm, often called the Louvain
 * algorithm, but it is faster and yields higher quality solutions. It can
 * optimize both modularity and the Constant Potts Model, which does not suffer
 * from the resolution-limit (see preprint http://arxiv.org/abs/1104.3083).
 *
 * </para><para>
 * The Leiden algorithm consists of three phases: (1) local moving of nodes,
 * (2) refinement of the partition and (3) aggregation of the network based on
 * the refined partition, using the non-refined partition to create an initial
 * partition for the aggregate network. In the local move procedure in the
 * Leiden algorithm, only nodes whose neighborhood has changed are visited. The
 * refinement is done by restarting from a singleton partition within each
 * cluster and gradually merging the subclusters. When aggregating, a single
 * cluster may then be represented by several nodes (which are the subclusters
 * identified in the refinement).
 *
 * </para><para>
 * The Leiden algorithm provides several guarantees. The Leiden algorithm is
 * typically iterated: the output of one iteration is used as the input for the
 * next iteration. At each iteration all clusters are guaranteed to be
 * connected and well-separated. After an iteration in which nothing has
 * changed, all nodes and some parts are guaranteed to be locally optimally
 * assigned. Finally, asymptotically, all subsets of all clusters are
 * guaranteed to be locally optimally assigned. For more details, please see
 * Traag, Waltman &amp; van Eck (2019).
 *
 * </para><para>
 * The objective function being optimized is
 *
 * </para><para>
 * 1 / 2m sum_ij (A_ij - gamma n_i n_j)d(s_i, s_j)
 *
 * </para><para>
 * where m is the total edge weight, A_ij is the weight of edge (i, j), gamma is
 * the so-called resolution parameter, n_i is the node weight of node i, s_i is
 * the cluster of node i and d(x, y) = 1 if and only if x = y and 0 otherwise.
 * By setting n_i = k_i, the degree of node i, and dividing gamma by 2m, you
 * effectively obtain an expression for modularity. Hence, the standard
 * modularity will be optimized when you supply the degrees as \c node_weights
 * and by supplying as a resolution parameter 1.0/(2*m), with m the number of
 * edges.
 *
 * \param graph The input graph. It must be an undirected graph.
 * \param edge_weights Numeric vector containing edge weights. If \c NULL, every edge
 *    has equal weight of 1. The weights need not be non-negative.
 * \param node_weights Numeric vector containing node weights.
 * \param resolution_parameter The resolution parameter used, which is
 *    represented by gamma in the objective function mentioned in the
 *    documentation.
 * \param beta The randomness used in the refinement step when merging. A small
 *    amount of randomness (\c beta = 0.01) typically works well.
 * \param start Start from membership vector. If this is true, the optimization
 *    will start from the provided membership vector. If this is false, the
 *    optimization will start from a singleton partition.
 * \param membership The membership vector. This is both used as the initial
 *    membership from which optimisation starts and is updated in place. It
 *    must hence be properly initialized. When finding clusters from scratch it
 *    is typically started using a singleton clustering. This can be achieved
 *    using \c igraph_vector_init_seq.
 * \param nb_clusters The number of clusters contained in \c membership. Must
 *    not be a \c NULL pointer.
 * \param quality The quality of the partition, in terms of the objective
 *    function as included in the documentation. If \c NULL the quality will
 *    not be calculated.
 * \return Error code.
 *
 * Time complexity: near linear on sparse graphs.
 *
 * \example examples/simple/igraph_community_leiden.c
 */
int igraph_community_leiden(const igraph_t *graph,
                            const igraph_vector_t *edge_weights, const igraph_vector_t *node_weights,
                            const igraph_real_t resolution_parameter, const igraph_real_t beta, const igraph_bool_t start,
                            igraph_vector_t *membership, igraph_integer_t *nb_clusters, igraph_real_t *quality) {
    igraph_vector_t *i_edge_weights, *i_node_weights;
    igraph_integer_t n = igraph_vcount(graph);

    if (start) {
        if (!membership) {
            IGRAPH_ERROR("Cannot start optimization if membership is missing", IGRAPH_EINVAL);
        }

        if (igraph_vector_size(membership) != n) {
            IGRAPH_ERROR("Initial membership length does not equal the number of vertices", IGRAPH_EINVAL);
        }
    } else {
        int i;
        if (!membership)
            IGRAPH_ERROR("Membership vector should be supplied and initialized, "
                         "even when not starting optimization from it", IGRAPH_EINVAL);

        igraph_vector_resize(membership, n);
        for (i = 0; i < n; i++) {
            VECTOR(*membership)[i] = i;
        }
    }


    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Leiden algorithm is only implemented for undirected graphs", IGRAPH_EINVAL);
    }

    /* Check edge weights to possibly use default */
    if (!edge_weights) {
        i_edge_weights = IGRAPH_CALLOC(1, igraph_vector_t);
        if (i_edge_weights == 0) {
            IGRAPH_ERROR("Leiden algorithm failed, could not allocate memory for edge weights", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, i_edge_weights);
        IGRAPH_CHECK(igraph_vector_init(i_edge_weights, igraph_ecount(graph)));
        IGRAPH_FINALLY(igraph_vector_destroy, i_edge_weights);
        igraph_vector_fill(i_edge_weights, 1);
    } else {
        i_edge_weights = (igraph_vector_t*)edge_weights;
    }

    /* Check edge weights to possibly use default */
    if (!node_weights) {
        i_node_weights = IGRAPH_CALLOC(1, igraph_vector_t);
        if (i_node_weights == 0) {
            IGRAPH_ERROR("Leiden algorithm failed, could not allocate memory for node weights", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, i_node_weights);
        IGRAPH_CHECK(igraph_vector_init(i_node_weights, n));
        IGRAPH_FINALLY(igraph_vector_destroy, i_node_weights);
        igraph_vector_fill(i_node_weights, 1);
    } else {
        i_node_weights = (igraph_vector_t*)node_weights;
    }

    /* Perform actual Leiden algorithm */
    IGRAPH_CHECK(igraph_i_community_leiden(graph, i_edge_weights, i_node_weights,
                                           resolution_parameter, beta,
                                           membership, nb_clusters, quality));

    if (!edge_weights) {
        igraph_vector_destroy(i_edge_weights);
        IGRAPH_FREE(i_edge_weights);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (!node_weights) {
        igraph_vector_destroy(i_node_weights);
        IGRAPH_FREE(i_node_weights);
        IGRAPH_FINALLY_CLEAN(2);
    }

    return IGRAPH_SUCCESS;
}
