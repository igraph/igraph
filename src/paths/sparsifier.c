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
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_random.h"

#include "core/math.h"

#include "config.h"

static void igraph_i_collect_lightest_edges_to_clusters(
    const igraph_adjlist_t *adjlist,
    const igraph_inclist_t *inclist,
    const igraph_vector_t *weights,
    const igraph_vector_t *clustering,
    const igraph_vector_bool_t *is_cluster_sampled,
    long int v, 
    igraph_vector_t *lightest_neighbor,
    igraph_vector_t *lightest_weight,
    long int *nearest_neighboring_sampled_cluster
) {
    // This internal function gets the residual graph, the clustering, the sampled clustering and
    // the vector and return the lightest edge to each neighboring cluster and the index of the lightest 
    // sampled cluster (if any)

    igraph_real_t lightest_weight_to_sampled = INFINITY;
    igraph_vector_int_t* adjacent_nodes = igraph_adjlist_get(adjlist, v);
    igraph_vector_int_t* incident_edges = igraph_inclist_get(inclist, v);
    long int i, nlen = igraph_vector_int_size(incident_edges);

    for (i = 0; i < nlen; i++) {
        long int neighbor_node = VECTOR(*adjacent_nodes)[i];
        long int edge = VECTOR(*incident_edges)[i];
        long int neighbor_cluster = VECTOR(*clustering)[neighbor_node];
        igraph_real_t weight = weights ? VECTOR(*weights)[edge] : 1;

        // If the weight of the edge being considered is smaller than the weight
        // of the lightest edge found so far that connects v to the same
        // cluster, remember the new minimum.
        if (VECTOR(*lightest_weight)[neighbor_cluster] > weight) {
            VECTOR(*lightest_weight)[neighbor_cluster] = weight;
            VECTOR(*lightest_neighbor)[neighbor_cluster] = edge;

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
}

/**
 * \ingroup structural
 * \function igraph_spanner
 * \brief Calculates a spanner of a graph with a given stretch factor.
 * 
 * A spanner of a graph G = (V,E) with a stretch t is a subgraph
 * H = (V,Es) such that Es is a subset of E and the distance 
 * between any pair of nodes in H is at most t times the distance
 * in G. The returned graph is always a spanner of the 
 * given graph with the specified stretch. For weighted graphs the
 * number of edges in the spanner is O(k * n^(1 + 1 / k)), where k is
 * k = (stretch + 1) / 2,  m is the number of edges and n is the number
 * of nodes in G. For unweighted graphs the number of edges is 
 * O(n^(1 + 1 / k) + kn).
 * 
 * </para><para>
 * This function is based on the algorithm of Baswana and Sen: "A Simple and
 * Linear Time Randomized Algorithm for Computing Sparse Spanners in 
 * Weighted Graphs"  
 *
 * \param graph An undirected connected graph object. If the graph
 *        is directed, the directions of the edges will be ignored.
 * \param spanenr An initialized vector, the IDs of the edges that constitute
 *        the calculated spanner will be returned here. Use
 *        \ref igraph_subgraph_edges() to extract the spanner as a separate
 *        graph object.
 * \param stretch The stretch factor of the spanner.
 * \param weights The edge weights or NULL. 
 * 
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \endclist
 *
 * Time complexity: The algorithm is a randomized Las Vegas algorithm. The expected
 *                  running time is O(km) where k is the value mentioned above.
 */
int igraph_spanner (const igraph_t *graph, igraph_vector_t *spanner,
        igraph_real_t stretch, igraph_vector_t *weights) { 

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int i, j, v, index, nlen, neighbor, cluster;
    double sample_prob, k = (stretch + 1) / 2, weight, lightest_sampled_weight;
    igraph_vector_t clustering, lightest_neighbor, lightest_weight;
    igraph_vector_bool_t is_cluster_sampled;
    igraph_vector_bool_t is_edge_in_spanner;
    igraph_vector_t new_clustering;
    igraph_vector_int_t *adjacent_vertices;
    igraph_vector_int_t *incident_edges;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    igraph_integer_t edge;

    if (spanner == 0) {
        return IGRAPH_SUCCESS;
    }

    /* Test validity of stretch factor */
    if (stretch < 1) {
        IGRAPH_ERROR("Stretch factor must be at least 1", IGRAPH_EINVAL);
    }

    /* Test validity of weights vector */
    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
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
    }

    // Clear the vector that will contain the IDs of the edges in the spanner
    igraph_vector_clear(spanner);

    // Create an incidence list reprsentation of the graph and also create the
    // corresponding adjacency list. The residual graph will not be constructed
    // explicitly; it will only exist in terms of the incidence and the adjacency
    // lists, maintained in parallel as the edges are removed from the residual
    // graph.
    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_OUT, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    IGRAPH_CHECK(igraph_adjlist_init_from_inclist(graph, &adjlist, &inclist));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    // Phase 1: forming the clusters
    // Create a vector which maps the nodes to the centers of the corresponding
    // clusters. At the beginning each node is its own cluster center.
    IGRAPH_CHECK(igraph_vector_init_seq(&clustering, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &clustering);

    // A mapping vector which indicates the neighboring edge with the smallest
    // weight for each cluster central
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_neighbor, no_of_nodes);

    // A mapping vector which indicated the minimum weight to each neighboring cluster
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_weight, no_of_nodes);

    IGRAPH_VECTOR_INIT_FINALLY(&new_clustering, no_of_nodes);

    // A boolean vector whose i-th element is 1 if the i-th vertex is a cluster
    // center that is sampled in the current iteration, 0 otherwise
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&is_cluster_sampled, no_of_nodes);

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&is_edge_in_spanner, no_of_edges);

    sample_prob = pow(no_of_nodes, -1 / k);

#define ADD_EDGE_TO_SPANNER {                                   \
    if (!VECTOR(is_edge_in_spanner)[edge]) {         \
        VECTOR(is_edge_in_spanner)[edge] = 1;        \
        IGRAPH_CHECK(igraph_vector_push_back(spanner, edge));   \
    }                                                           \
}

    for (i = 0; i < k - 1; i++) {
        igraph_vector_fill(&new_clustering, -1);
        igraph_vector_bool_fill(&is_cluster_sampled, 0);

        // Step 1: sample cluster centers
        RNG_BEGIN();
        for (j = 0; j < no_of_nodes; j++) {
            if (VECTOR(clustering)[j] == j && RNG_UNIF01() < sample_prob) {
                VECTOR(is_cluster_sampled)[j] = 1;
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
            
            igraph_vector_fill(&lightest_neighbor, -1);
            igraph_vector_fill(&lightest_weight, INFINITY);
            
            // Step 2: find the lightest edge that connects vertex v to its
            // neighboring sampled clusters
            long int nearest_neighboring_sampled_cluster = -1;
            igraph_i_collect_lightest_edges_to_clusters(
                &adjlist,
                &inclist,
                weights,
                &clustering,
                &is_cluster_sampled,
                v,
                &lightest_neighbor,
                &lightest_weight,
                &nearest_neighboring_sampled_cluster
            );

            // Step 3: add edges to spanner
            if (nearest_neighboring_sampled_cluster == -1) {
                // Case 3(a) from the paper: v is not adjacent to any of the
                // sampled clusters.

                // Add lightest edge which connects vertex v to each neighboring
                // cluster (none of which are sampled) 
                for (j = 0; j < no_of_nodes; j++) {
                    edge = VECTOR(lightest_neighbor)[j];
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
                edge = VECTOR(lightest_neighbor)[nearest_neighboring_sampled_cluster];
                ADD_EDGE_TO_SPANNER;

                // 'lightest_sampled_weight' is the weight of the lightest edge connecting v to
                // one of the sampled clusters. This is where v will belong in
                // the new clustering.
                lightest_sampled_weight = VECTOR(lightest_weight)[nearest_neighboring_sampled_cluster];
                VECTOR(new_clustering)[v] = nearest_neighboring_sampled_cluster;

                // Add to the spanner light edges with weight less than 'lightest_sampled_weight'
                for (j = 0; j < no_of_nodes; j++) {
                    if (VECTOR(lightest_weight)[j] < lightest_sampled_weight) {
                        edge = VECTOR(lightest_neighbor)[j];
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
        }

        // Commit the new clustering
        igraph_vector_update(&clustering, &new_clustering);

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
            igraph_vector_fill(&lightest_neighbor, -1);
            igraph_vector_fill(&lightest_weight, INFINITY);
            igraph_i_collect_lightest_edges_to_clusters(
                &adjlist,
                &inclist,
                weights,
                &clustering,
                /* is_cluster_sampled = */ NULL,
                v, 
                &lightest_neighbor,
                &lightest_weight,
                NULL
            );
            for (j = 0; j < no_of_nodes; j++) {
                edge = VECTOR(lightest_neighbor)[j];
                if (edge != -1) {
                    ADD_EDGE_TO_SPANNER;
                }
            }
        }
    }

    // Free memory
    igraph_vector_bool_destroy(&is_edge_in_spanner);
    igraph_inclist_destroy(&inclist);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_destroy(&clustering);
    igraph_vector_destroy(&lightest_neighbor);
    igraph_vector_destroy(&lightest_weight);
    igraph_vector_destroy(&new_clustering);
    igraph_vector_bool_destroy(&is_cluster_sampled);
    IGRAPH_FINALLY_CLEAN(8);

    return IGRAPH_SUCCESS;
}
