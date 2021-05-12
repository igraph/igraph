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
#include "igraph_heap.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_random.h"

#include "core/math.h"

#include "config.h"

static void igraph_i_collect_lightest_edges_to_clusters(
    const igraph_t *residual_graph,
    const igraph_vector_t *weights,
    const igraph_vector_t *clustering,
    const igraph_vector_bool_t *is_cluster_sampled,
    long int v, 
    const igraph_vector_int_t *neis,
    igraph_vector_t *lightest_neighbor,
    igraph_vector_t *lightest_weight,
    long int *nearest_neighboring_sampled_cluster
) {
    // This internal function gets the residual graph, the clustering, the sampled clustering and
    // the vector and return the lightest edge to each neighboring cluster and the index of the lightest 
    // sampled cluster (if any)

    float lightest_weight_to_sampled = INFINITY;
    long int i, nlen = igraph_vector_int_size(neis);

    for (i = 0; i < nlen; i++) {
        long int edge = VECTOR(*neis)[i];
        long int neighbor_node = IGRAPH_OTHER(residual_graph, edge, v);
        long int neighbor_cluster = VECTOR(*clustering)[neighbor_node];
        float weight = VECTOR(*weights)[edge];

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
    long int i = 0, nlen, neighbor, cluster_of_v;
    long int edge, edge_old;
    double sample_prob, size_limit, k = (stretch + 1) / 2;
    igraph_t residual_graph;
    igraph_vector_t clustering, centers, lightest_neighbor, lightest_weight, residual_weight;
    igraph_vector_bool_t is_cluster_sampled;
    igraph_vector_t new_clustering, eids;
    igraph_vector_int_t *neis;
    igraph_integer_t from, to;
    igraph_integer_t edge_in_original_graph;
    igraph_inclist_t inclist;
    igraph_es_t es;
    igraph_heap_t edges_to_remove, edges_to_add;

    if (spanner == 0) {
        return IGRAPH_SUCCESS;
    }

    /* Test validity of stretch factor */
    if (stretch < 1) {
        IGRAPH_ERROR("Stretch factor must be at least 1", IGRAPH_EINVAL);
    }

    // If the weights is empty (the graph is unweighted) then assign a weight vector of ones
    if (!weights) {
        IGRAPH_VECTOR_INIT_FINALLY(&residual_weight, no_of_edges);
        igraph_vector_fill(&residual_weight, 1);
    } else {
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
        IGRAPH_CHECK(igraph_vector_copy(&residual_weight, weights));
        IGRAPH_FINALLY(igraph_vector_destroy, &residual_weight);
    }

    // Clear the vector that will contain the IDs of the edges in the spanner
    igraph_vector_clear(spanner);

    // Create a residual graph with V' = V and E' = E. If the graph is unweighted,
    // then each edge has weight of 1.
    IGRAPH_CHECK(igraph_copy(&residual_graph, graph));
    IGRAPH_FINALLY(igraph_destroy, &residual_graph);

    // Phase 1: forming the clusters
    // Create a vector which maps the nodes to the centers of the corresponding
    // clusters. At the beginning each node is its own cluster center.
    IGRAPH_CHECK(igraph_vector_init_seq(&clustering, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &clustering);

    // A vector with the list of the clustering centers. In the beginning of the 
    // algorithm all the nodes are centers to their clusters. This is simply
    // an alternative representation of 'clustering' but we need it to be able
    // to sample efficiently from the current set of centers.
    IGRAPH_CHECK(igraph_vector_init_seq(&centers, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &centers);

    // A mapping vector which indicates the neighboring edge with the smallest
    // weight for each cluster central
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_neighbor, no_of_nodes);

    // A mapping vector which indicated the minimum weight to each neighboring cluster
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_weight, no_of_nodes);

    IGRAPH_VECTOR_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&new_clustering, no_of_nodes);

    // A boolean vector whose i-th element is 1 if the i-th vertex is a cluster
    // center that is sampled in the current iteration, 0 otherwise
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&is_cluster_sampled, no_of_nodes);

    IGRAPH_CHECK(igraph_heap_init(&edges_to_add, 0));
    IGRAPH_FINALLY(igraph_heap_destroy, &edges_to_add);
    IGRAPH_CHECK(igraph_heap_reserve(&edges_to_add, no_of_edges));

    IGRAPH_CHECK(igraph_heap_init(&edges_to_remove, 0));
    IGRAPH_FINALLY(igraph_heap_destroy, &edges_to_remove);
    IGRAPH_CHECK(igraph_heap_reserve(&edges_to_remove, no_of_edges));

    sample_prob = pow(no_of_nodes, -1 / k);
    size_limit = pow(no_of_nodes, 1+1/k);

    while (i < k-1) {
        IGRAPH_CHECK(igraph_inclist_init(&residual_graph, &inclist, IGRAPH_OUT, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

        igraph_vector_fill(&new_clustering, -1);
        igraph_vector_bool_fill(&is_cluster_sampled, 0);

        // Step 1: sample centers
        RNG_BEGIN();
        int size_centers = igraph_vector_size(&centers);
        for (int j = 0; j < size_centers; j++) {
            if (RNG_UNIF01() < sample_prob) {
                long int center = VECTOR(centers)[j];
                VECTOR(is_cluster_sampled)[center] = 1;
            }
        }
        RNG_END();

        // Step 2 and 3
        for (long int v = 0; v < no_of_nodes; v++) {
            // If v is inside a cluster and the cluster of v is sampled, then continue
            cluster_of_v = VECTOR(clustering)[v];
            if (cluster_of_v != -1 && VECTOR(is_cluster_sampled)[cluster_of_v]) {
                VECTOR(new_clustering)[v] = cluster_of_v;
                continue;
            }
            
            neis = igraph_inclist_get(&inclist, v);
            nlen = igraph_vector_int_size(neis);
            igraph_vector_fill(&lightest_neighbor, -1);
            igraph_vector_fill(&lightest_weight, INFINITY);
            
            // Step 2: find the lightest edge that connects vertex v to its
            // neighboring sampled clusters
            long int nearest_neighboring_sampled_cluster = -1;
            igraph_i_collect_lightest_edges_to_clusters(
                &residual_graph, 
                &residual_weight,
                &clustering,
                &is_cluster_sampled,
                v, 
                neis, 
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
                for (long int j = 0; j < no_of_nodes; j++) {
                    if (VECTOR(lightest_neighbor)[j] != -1) {
                        edge = VECTOR(lightest_neighbor)[j];
                        IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                    }
                }

                // Remove all edges incident on v from the graph
                for (int j = 0; j < nlen; j++) {
                    edge = VECTOR(*neis)[j];
                    IGRAPH_CHECK(igraph_heap_push(&edges_to_remove, (igraph_real_t) edge));
                }
            } else {
                // Case 3(b) from the paper: v is adjacent to at least one of
                // the sampled clusters

                // add the edge connecting to the lightest sampled cluster
                edge = VECTOR(lightest_neighbor)[nearest_neighboring_sampled_cluster];
                IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                double lightest_sample_weight = VECTOR(lightest_weight)[nearest_neighboring_sampled_cluster];
                VECTOR(new_clustering)[v] = nearest_neighboring_sampled_cluster;

                // Add to the spanner light edges with weight less then the lightest_sample_weight
                for (int j = 0; j < no_of_nodes; j++) {
                    if (VECTOR(lightest_weight)[j] < lightest_sample_weight) {
                        edge = VECTOR(lightest_neighbor)[j];
                        IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                    }
                }

                // Remove edges to centers with edge weight less than lightest_sample_weight
                for (int j = 0; j < nlen; j++) {
                    neighbor = IGRAPH_OTHER(&residual_graph, VECTOR(*neis)[j], v);
                    long int neighbor_cluster = VECTOR(clustering)[neighbor];
                    double weight = VECTOR(lightest_weight)[neighbor_cluster];
                    if ((neighbor_cluster == nearest_neighboring_sampled_cluster) || (weight < lightest_sample_weight)) {
                        edge = VECTOR(*neis)[j];
                        IGRAPH_CHECK(igraph_heap_push(&edges_to_remove, (igraph_real_t) edge));
                    }
                }
            }
        }

        // Check if this iteration added too many edges to the spanner
        if (igraph_heap_size(&edges_to_add) > size_limit) {
            igraph_inclist_destroy(&inclist);
            IGRAPH_FINALLY_CLEAN(1);
            igraph_heap_clear(&edges_to_add);
            igraph_heap_clear(&edges_to_remove);
            continue;
        }
        i = i+1;

        // Add the edges to the spanner
        edge_old = -1;
        while (!igraph_heap_empty(&edges_to_add)) {
            edge = (long int) igraph_heap_delete_top(&edges_to_add);
            if (edge != edge_old) {
                // only if the edge is new, meaning it wasn't inserted by the j-1 element in the list
                // then insert the edge to the spanner. 
                igraph_edge(&residual_graph, edge, &from, &to);
                IGRAPH_CHECK(
                    igraph_get_eid(graph, &edge_in_original_graph, from, to, /* directed = */ 0, /* error = */ 1)
                );
                IGRAPH_CHECK(igraph_vector_push_back(spanner, edge_in_original_graph));
            }
            edge_old = edge;
        }

        // Delete edges from residual graph
        edge_old = -1;
        while (!igraph_heap_empty(&edges_to_remove)) {
            // Remove the weights in the residual_weight vector in nondecreasing order
            edge = (long int) igraph_heap_delete_top(&edges_to_remove);
            if (edge != edge_old) {
                igraph_es_1(&es, edge);
                IGRAPH_CHECK(igraph_delete_edges(&residual_graph, es));
                igraph_vector_remove_section(&residual_weight, edge, edge+1);
                igraph_es_destroy(&es);
            }
            edge_old = edge;
        }

        // Copy old clustering to new clutering
        igraph_vector_fill(&clustering, -1);
        igraph_vector_resize(&centers, 0);
        for (int v = 0; v < no_of_nodes; v++) {
            long int new_cluster = VECTOR(new_clustering)[v];
            if (new_cluster == v) {
                igraph_vector_push_back(&centers, new_cluster); // Mark that vertex new_cluster is a center
            }
            VECTOR(clustering)[v] = new_cluster;
        }

        // Remove intra-cluster edges
        for (int v = 0; v < no_of_nodes; v++) {
            igraph_vector_clear(&eids);
            IGRAPH_CHECK(igraph_incident(&residual_graph, &eids, v, IGRAPH_OUT));
            nlen = igraph_vector_size(&eids);
            for (int j = nlen-1; j >= 0; j--) {
                edge = VECTOR(eids)[j];
                long int neighbor = IGRAPH_OTHER(&residual_graph, edge, v);
                if (VECTOR(clustering)[neighbor] == VECTOR(clustering)[v]) {
                    igraph_es_1(&es, edge);
                    IGRAPH_CHECK(igraph_delete_edges(&residual_graph, es));
                    igraph_vector_remove_section(&residual_weight, edge, edge+1);
                    igraph_es_destroy(&es);
                }
            }
        }

        // Free
        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    // Phase 2: vertex_clustering joining
    IGRAPH_CHECK(igraph_inclist_init(&residual_graph, &inclist, IGRAPH_OUT, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    for (int v = 0; v < no_of_nodes; v++) {
        if (VECTOR(clustering)[v] != -1) {
            igraph_vector_fill(&lightest_neighbor, -1);
            igraph_vector_fill(&lightest_weight, INFINITY);
            neis = igraph_inclist_get(&inclist, v);
            igraph_i_collect_lightest_edges_to_clusters(
                &residual_graph, 
                &residual_weight,
                &clustering,
                /* is_cluster_sampled = */ NULL,
                v, 
                neis,
                &lightest_neighbor,
                &lightest_weight,
                NULL
            );
            for (int j = 0; j < no_of_nodes; j++) {
                if (VECTOR(lightest_neighbor)[j] != -1) {
                    edge = VECTOR(lightest_neighbor)[j];
                    IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                }
            }
        }
    }

    edge_old = -1;
    while (!igraph_heap_empty(&edges_to_add)) {
        edge = (long int)igraph_heap_delete_top(&edges_to_add);
        if (edge != edge_old) {
            // only if the edge is new, meaning it wasn't inserted by the j-1 element in the list
            // then insert the edge to the spanner. 
            igraph_edge(&residual_graph, edge, &from, &to);
            IGRAPH_CHECK(
                igraph_get_eid(graph, &edge_in_original_graph, from, to, /* directed = */ 0, /* error = */ 1)
            );
            IGRAPH_CHECK(igraph_vector_push_back(spanner, edge_in_original_graph));
        }
        edge_old = edge;
    }
    
    // Free memory
    igraph_vector_destroy(&residual_weight);
    igraph_destroy(&residual_graph);
    igraph_vector_destroy(&clustering);
    igraph_vector_destroy(&centers);
    igraph_inclist_destroy(&inclist);
    igraph_vector_destroy(&lightest_neighbor);
    igraph_vector_destroy(&lightest_weight);
    igraph_vector_destroy(&new_clustering);
    igraph_vector_bool_destroy(&is_cluster_sampled);
    igraph_heap_destroy(&edges_to_remove);
    igraph_heap_destroy(&edges_to_add);
    igraph_vector_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(12);

    return IGRAPH_SUCCESS;
}
