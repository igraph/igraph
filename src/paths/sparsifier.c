/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_random.h"

#include "core/interruption.h"

/*
 * This internal function gets the adjacency and incidence list representation
 * of the current residual graph, the weight vector, the current assignment of
 * the vertices to clusters, whether the i-th cluster is sampled, and the
 * index of a single node v. The function updates the given lightest_eid vector
 * such that the i-th element contains the ID of the lightest edge that leads
 * from node v to cluster i. Similarly, the lightest_weight vector is updated
 * to contain the weights of these edges.
 *
 * When the is_cluster_sampled vector is provided, the
 * nearest_neighboring_sampled_cluster pointer is also updated to the index of
 * the cluster that has the smallest weight among the _sampled_ ones.
 *
 * As a pre-condition, this function requires the lightest_eid vector to be
 * filled with -1 and the lightest_weight vector to be filled with infinity.
 * This is _not_ checked within the function.
 *
 * Use the igraph_i_clean_lightest_edges_to_clusters() function to clear these vectors
 * after you are done with them. Avoid using igraph_vector_fill() because that
 * one is O(|V|), while igraph_i_clean_lightest_edge_vector() is O(d) where d
 * is the degree of the vertex.
 */
static igraph_error_t igraph_i_collect_lightest_edges_to_clusters(
    const igraph_adjlist_t *adjlist,
    const igraph_inclist_t *inclist,
    const igraph_vector_t *weights,
    const igraph_vector_int_t *clustering,
    const igraph_vector_bool_t *is_cluster_sampled,
    igraph_integer_t v,
    igraph_vector_int_t *lightest_eid,
    igraph_vector_t *lightest_weight,
    igraph_vector_int_t *dirty_vids,
    igraph_integer_t *nearest_neighboring_sampled_cluster
) {
    // This internal function gets the residual graph, the clustering, the sampled clustering and
    // the vector and return the lightest edge to each neighboring cluster and the index of the lightest
    // sampled cluster (if any)

    igraph_real_t lightest_weight_to_sampled = INFINITY;
    igraph_vector_int_t* adjacent_nodes = igraph_adjlist_get(adjlist, v);
    igraph_vector_int_t* incident_edges = igraph_inclist_get(inclist, v);
    igraph_integer_t i, nlen = igraph_vector_int_size(incident_edges);

    for (i = 0; i < nlen; i++) {
        igraph_integer_t neighbor_node = VECTOR(*adjacent_nodes)[i];
        igraph_integer_t edge = VECTOR(*incident_edges)[i];
        igraph_integer_t neighbor_cluster = VECTOR(*clustering)[neighbor_node];
        igraph_real_t weight = weights ? VECTOR(*weights)[edge] : 1;

        // If the weight of the edge being considered is smaller than the weight
        // of the lightest edge found so far that connects v to the same
        // cluster, remember the new minimum.
        if (VECTOR(*lightest_weight)[neighbor_cluster] > weight) {
            VECTOR(*lightest_weight)[neighbor_cluster] = weight;
            VECTOR(*lightest_eid)[neighbor_cluster] = edge;

            IGRAPH_CHECK(igraph_vector_int_push_back(dirty_vids, neighbor_cluster));

            // Also, if this cluster happens to be a sampled cluster, also update
            // the variables that store which is the lightest edge that connects
            // v to any of the sampled clusters.
            if (is_cluster_sampled) {
                if ((VECTOR(*is_cluster_sampled)[neighbor_cluster]) && (lightest_weight_to_sampled > weight)) {
                    lightest_weight_to_sampled = weight;
                    *nearest_neighboring_sampled_cluster = neighbor_cluster;
                }
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static void igraph_i_clear_lightest_edges_to_clusters(
    igraph_vector_int_t *dirty_vids,
    igraph_vector_int_t *lightest_eid,
    igraph_vector_t *lightest_weight
) {
    igraph_integer_t i, n = igraph_vector_int_size(dirty_vids);
    for (i = 0; i < n; i++) {
        igraph_integer_t vid = VECTOR(*dirty_vids)[i];
        VECTOR(*lightest_weight)[vid] = INFINITY;
        VECTOR(*lightest_eid)[vid] = -1;
    }

    igraph_vector_int_clear(dirty_vids);
}

/**
 * \ingroup structural
 * \function igraph_spanner
 * \brief Calculates a spanner of a graph with a given stretch factor.
 *
 * A spanner of a graph <code>G = (V,E)</code> with a stretch \c t is a
 * subgraph <code>H = (V,Es)</code> such that \c Es is a subset of \c E
 * and the distance between any pair of nodes in \c H is at most \c t
 * times the distance in \c G. The returned graph is always a spanner of
 * the given graph with the specified stretch. For weighted graphs the
 * number of edges in the spanner is <code>O(k n^(1 + 1 / k))</code>, where
 * \c k is <code>k = (t + 1) / 2</code>,  \c m is the number of edges
 * and \c n is the number of nodes in \c G. For unweighted graphs the number
 * of edges is <code>O(n^(1 + 1 / k) + kn)</code>.
 *
 * </para><para>
 * This function is based on the algorithm of Baswana and Sen: "A Simple and
 * Linear Time Randomized Algorithm for Computing Sparse Spanners in
 * Weighted Graphs". https://doi.org/10.1002/rsa.20130
 *
 * \param graph An undirected connected graph object. If the graph
 *        is directed, the directions of the edges will be ignored.
 * \param spanner An initialized vector, the IDs of the edges that constitute
 *        the calculated spanner will be returned here. Use
 *        \ref igraph_subgraph_from_edges() to extract the spanner as a separate
 *        graph object.
 * \param stretch The stretch factor \c t of the spanner.
 * \param weights The edge weights or \c NULL.
 *
 * \return Error code:
 * \clist
 * \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 * \endclist
 *
 * Time complexity: The algorithm is a randomized Las Vegas algorithm. The expected
 * running time is O(km) where k is the value mentioned above and m is the number
 * of edges.
 */
igraph_error_t igraph_spanner(const igraph_t *graph, igraph_vector_int_t *spanner,
        igraph_real_t stretch, const igraph_vector_t *weights) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t i, j, v, nlen, neighbor, cluster;
    igraph_real_t sample_prob, k = (stretch + 1) / 2, weight, lightest_sampled_weight;
    igraph_vector_int_t clustering, lightest_eid;
    igraph_vector_t lightest_weight;
    igraph_vector_bool_t is_cluster_sampled;
    igraph_vector_bool_t is_edge_in_spanner;
    igraph_vector_int_t new_clustering;
    igraph_vector_int_t dirty_vids;
    igraph_vector_int_t *adjacent_vertices;
    igraph_vector_int_t *incident_edges;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    igraph_integer_t edge;
    igraph_integer_t index;

    if (spanner == NULL) {
        return IGRAPH_SUCCESS;
    }

    /* Test validity of stretch factor */
    if (stretch < 1) {
        IGRAPH_ERROR("Stretch factor must be at least 1.", IGRAPH_EINVAL);
    }

    /* Test validity of weights vector */
    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match.", IGRAPH_EINVAL);
        }
        if (no_of_edges > 0) {
            igraph_real_t min = igraph_vector_min(weights);
            if (min < 0) {
                IGRAPH_ERROR("Weight vector must be non-negative.", IGRAPH_EINVAL);
            }
            else if (isnan(min)) {
                IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
            }
        }
    }

    // Clear the vector that will contain the IDs of the edges in the spanner
    igraph_vector_int_clear(spanner);

    // Create an incidence list representation of the graph and also create the
    // corresponding adjacency list. The residual graph will not be constructed
    // explicitly; it will only exist in terms of the incidence and the adjacency
    // lists, maintained in parallel as the edges are removed from the residual
    // graph.
    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_adjlist_init_from_inclist(graph, &adjlist, &inclist));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    // Phase 1: forming the clusters
    // Create a vector which maps the nodes to the centers of the corresponding
    // clusters. At the beginning each node is its own cluster center.
    IGRAPH_CHECK(igraph_vector_int_init_range(&clustering, 0, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &clustering);

    // A mapping vector which indicates the neighboring edge with the smallest
    // weight for each cluster central, for a single vertex of interest.
    // Preconditions needed by igraph_i_collect_lightest_edges_to_clusters()
    // are enforced here.
    IGRAPH_VECTOR_INT_INIT_FINALLY(&lightest_eid, no_of_nodes);
    igraph_vector_int_fill(&lightest_eid, -1);

    // A mapping vector which indicated the minimum weight to each neighboring
    // cluster, for a single vertex of interest.
    // Preconditions needed by igraph_i_collect_lightest_edges_to_clusters()
    // are enforced here.
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_weight, no_of_nodes);
    igraph_vector_fill(&lightest_weight, IGRAPH_INFINITY);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&new_clustering, no_of_nodes);

    // A boolean vector whose i-th element is 1 if the i-th vertex is a cluster
    // center that is sampled in the current iteration, 0 otherwise
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&is_cluster_sampled, no_of_nodes);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&is_edge_in_spanner, no_of_edges);

    // Temporary vector used by igraph_i_collect_lightest_edges_to_clusters()
    // to keep track of the nodes that it has written to
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dirty_vids, 0);

    sample_prob = pow(no_of_nodes, -1 / k);

#define ADD_EDGE_TO_SPANNER \
    if (!VECTOR(is_edge_in_spanner)[edge]) { \
        VECTOR(is_edge_in_spanner)[edge] = true; \
        IGRAPH_CHECK(igraph_vector_int_push_back(spanner, edge)); \
    }

    igraph_vector_fill(&lightest_weight, INFINITY);

    for (i = 0; i < k - 1; i++) {
        IGRAPH_ALLOW_INTERRUPTION();

        igraph_vector_int_fill(&new_clustering, -1);
        igraph_vector_bool_fill(&is_cluster_sampled, false);

        // Step 1: sample cluster centers
        RNG_BEGIN();
        for (j = 0; j < no_of_nodes; j++) {
            if (VECTOR(clustering)[j] == j && RNG_UNIF01() < sample_prob) {
                VECTOR(is_cluster_sampled)[j] = true;
            }
        }
        RNG_END();

        // Step 2 and 3
        for (v = 0; v < no_of_nodes; v++) {
            // If v is inside a cluster and the cluster of v is sampled, then continue
            cluster = VECTOR(clustering)[v];
            if (cluster != -1 && VECTOR(is_cluster_sampled)[cluster]) {
                VECTOR(new_clustering)[v] = cluster;
                continue;
            }

            // Step 2: find the lightest edge that connects vertex v to its
            // neighboring sampled clusters
            igraph_integer_t nearest_neighboring_sampled_cluster = -1;
            IGRAPH_CHECK(igraph_i_collect_lightest_edges_to_clusters(
                &adjlist,
                &inclist,
                weights,
                &clustering,
                &is_cluster_sampled,
                v,
                &lightest_eid,
                &lightest_weight,
                &dirty_vids,
                &nearest_neighboring_sampled_cluster
            ));

            // Step 3: add edges to spanner
            if (nearest_neighboring_sampled_cluster == -1) {
                // Case 3(a) from the paper: v is not adjacent to any of the
                // sampled clusters.

                // Add lightest edge which connects vertex v to each neighboring
                // cluster (none of which are sampled)
                for (j = 0; j < no_of_nodes; j++) {
                    edge = VECTOR(lightest_eid)[j];
                    if (edge != -1) {
                        ADD_EDGE_TO_SPANNER;
                    }
                }

                // Remove all edges incident on v from the graph. Note that each
                // edge being removed occurs twice in the adjacency / incidence
                // lists
                adjacent_vertices = igraph_adjlist_get(&adjlist, v);
                incident_edges = igraph_inclist_get(&inclist, v);
                nlen = igraph_vector_int_size(incident_edges);
                for (j = 0; j < nlen; j++) {
                    neighbor = VECTOR(*adjacent_vertices)[j];
                    if (neighbor == v) {
                        /* should not happen as we did not ask for loop edges in
                         * the adjacency / incidence lists, but let's be defensive */
                        continue;
                    }

                    if (igraph_vector_int_search(
                        igraph_inclist_get(&inclist, neighbor),
                        0,
                        VECTOR(*incident_edges)[j],
                        &index
                    )) {
                        igraph_vector_int_remove_fast(igraph_adjlist_get(&adjlist, neighbor), index);
                        igraph_vector_int_remove_fast(igraph_inclist_get(&inclist, neighbor), index);
                    }
                }
                igraph_vector_int_clear(adjacent_vertices);
                igraph_vector_int_clear(incident_edges);
            } else {
                // Case 3(b) from the paper: v is adjacent to at least one of
                // the sampled clusters

                // add the edge connecting to the lightest sampled cluster
                edge = VECTOR(lightest_eid)[nearest_neighboring_sampled_cluster];
                ADD_EDGE_TO_SPANNER;

                // 'lightest_sampled_weight' is the weight of the lightest edge connecting v to
                // one of the sampled clusters. This is where v will belong in
                // the new clustering.
                lightest_sampled_weight = VECTOR(lightest_weight)[nearest_neighboring_sampled_cluster];
                VECTOR(new_clustering)[v] = nearest_neighboring_sampled_cluster;

                // Add to the spanner light edges with weight less than 'lightest_sampled_weight'
                for (j = 0; j < no_of_nodes; j++) {
                    if (VECTOR(lightest_weight)[j] < lightest_sampled_weight) {
                        edge = VECTOR(lightest_eid)[j];
                        ADD_EDGE_TO_SPANNER;
                    }
                }

                // Remove edges to centers with edge weight less than 'lightest_sampled_weight'
                adjacent_vertices = igraph_adjlist_get(&adjlist, v);
                incident_edges = igraph_inclist_get(&inclist, v);
                nlen = igraph_vector_int_size(incident_edges);
                for (j = 0; j < nlen; j++) {
                    neighbor = VECTOR(*adjacent_vertices)[j];
                    if (neighbor == v) {
                        /* should not happen as we did not ask for loop edges in
                         * the adjacency / incidence lists, but let's be defensive */
                        continue;
                    }

                    cluster = VECTOR(clustering)[neighbor];
                    weight = VECTOR(lightest_weight)[cluster];
                    if ((cluster == nearest_neighboring_sampled_cluster) || (weight < lightest_sampled_weight)) {
                        edge = VECTOR(*incident_edges)[j];

                        if (igraph_vector_int_search(
                            igraph_inclist_get(&inclist, neighbor), 0, edge, &index
                        )) {
                            igraph_vector_int_remove_fast(igraph_adjlist_get(&adjlist, neighbor), index);
                            igraph_vector_int_remove_fast(igraph_inclist_get(&inclist, neighbor), index);
                        }

                        igraph_vector_int_remove_fast(adjacent_vertices, j);
                        igraph_vector_int_remove_fast(incident_edges, j);

                        j--;
                        nlen--;
                    }
                }
            }

            // We don't need lightest_eids and lightest_weights any more so
            // clear them in O(d) time
            igraph_i_clear_lightest_edges_to_clusters(
                &dirty_vids, &lightest_eid, &lightest_weight
            );
        }

        // Commit the new clustering
        igraph_vector_int_update(&clustering, &new_clustering); /* reserved */

        // Remove intra-cluster edges
        for (v = 0; v < no_of_nodes; v++) {
            adjacent_vertices = igraph_adjlist_get(&adjlist, v);
            incident_edges = igraph_inclist_get(&inclist, v);
            nlen = igraph_vector_int_size(incident_edges);
            for (j = 0; j < nlen; j++) {
                neighbor = VECTOR(*adjacent_vertices)[j];
                edge = VECTOR(*incident_edges)[j];

                if (VECTOR(clustering)[neighbor] == VECTOR(clustering)[v]) {
                    /* We don't need to bother with removing the other copy
                     * of the edge from the incidence lists (and the corresponding
                     * vertices from the adjacency lists) because we will find
                     * them anyway as we are iterating over all nodes */
                    igraph_vector_int_remove_fast(adjacent_vertices, j);
                    igraph_vector_int_remove_fast(incident_edges, j);
                    j--;
                    nlen--;
                }
            }
        }
    }

    // Phase 2: vertex_clustering joining
    for (v = 0; v < no_of_nodes; v++) {
        if (VECTOR(clustering)[v] != -1) {
            IGRAPH_CHECK(igraph_i_collect_lightest_edges_to_clusters(
                &adjlist,
                &inclist,
                weights,
                &clustering,
                /* is_cluster_sampled = */ NULL,
                v,
                &lightest_eid,
                &lightest_weight,
                &dirty_vids,
                NULL
            ));
            for (j = 0; j < no_of_nodes; j++) {
                edge = VECTOR(lightest_eid)[j];
                if (edge != -1) {
                    ADD_EDGE_TO_SPANNER;
                }
            }
            igraph_i_clear_lightest_edges_to_clusters(&dirty_vids, &lightest_eid, &lightest_weight);
        }
    }

    // Free memory
    igraph_vector_int_destroy(&dirty_vids);
    igraph_vector_bool_destroy(&is_edge_in_spanner);
    igraph_vector_bool_destroy(&is_cluster_sampled);
    igraph_vector_int_destroy(&new_clustering);
    igraph_vector_destroy(&lightest_weight);
    igraph_vector_int_destroy(&lightest_eid);
    igraph_vector_int_destroy(&clustering);
    igraph_adjlist_destroy(&adjlist);
    igraph_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}
