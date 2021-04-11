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

static int igraph_i_lightest_edge_cluster(const igraph_t *residual_graph,
                                    const igraph_inclist_t *inclist, 
                                    long int v, 
                                    long int *lightest_sample_clustering_index,
                                    const igraph_vector_t *weights,
                                    const igraph_vector_t *clustering,
                                    const igraph_vector_t *sample_clustering,
                                    igraph_vector_t *lightest_neighbor,
                                    igraph_vector_t *lightest_weight) {
    // This internal function gets the residual graph, the clustering, the sampled clustering and
    // the vector and return the lightest edge to each neighboring cluster and the index of the lightest 
    // sampled cluster (if any)
    float lightest_weight_to_sampled = INFINITY;
    igraph_vector_int_t *neis = igraph_inclist_get(inclist, v);
    long int nlen = igraph_vector_int_size(neis);
    for (int i = 0; i < nlen; i++) {
        long int edge = VECTOR(*neis)[i];
        long int neighbor_node = IGRAPH_OTHER(residual_graph, edge, v);
        long int neighbor_cluster = VECTOR(*clustering)[neighbor_node];
        float weight = VECTOR(*weights)[edge];
        if (VECTOR(*lightest_weight)[neighbor_cluster] > weight) {
            VECTOR(*lightest_weight)[neighbor_cluster] = weight;
            VECTOR(*lightest_neighbor)[neighbor_cluster] = edge;

            // Find the lightest sampled cluster 
            if (sample_clustering) {
                if ((VECTOR(*sample_clustering)[neighbor_cluster]) && (lightest_weight_to_sampled > weight)) {
                    lightest_weight_to_sampled = weight;
                    *lightest_sample_clustering_index = neighbor_cluster;
                }
            }
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_spanner
 * \brief Retrun a spanner of a graph with a given stretch factor.
 * 
 * A spanner of a graph G = (V,E) with a stretch t is a subgraph
 * H = (V,Es) such that Es is a subset of E and the distance 
 * between any pair of nodes in H is at most t times the distance
 * is G. The returned graph is always a spanner of the 
 * given graph with the specified stretch. For weighted graphs the
 * number of edges in the spanner is O(k * n^(1 + 1 / k)) where k is
 * k = (stretch + 1) / 2,  m is the number of edges and n is the number
 * of nodes in G. For unweighted graphs the number of edges is 
 * O(n^(1 + 1 / k) + kn).
 * 
 * </para><para>
 * This function is based on Baswana and Sen random algorithm: "A Simple and
 * Linear Time Randomized Algorithm for Computing Sparse Spanners in 
 * Weighted Graphs"  
 *
 * \param graph An undirected connected graph object. If the graph
 *        is directed, the directions of the edges would be ignored.
 * \param spanner Pointer to an uninitalized graph_t pointer. The
 *        function would return the spanner in this pointer
 * \param stretch The stretch factor of the spanner.
 * \param weights The edge weights or NULL. 
 * \param spanner_weight Pointer to uninitalized vector or NULL. If not NULL,
 *        the function would save the weights of the edges in the spanner
 *        accourding to the edges IDs in the spanner graph.
 * \param seed A pointer to unsigned long int or NULL. If not NULL would be
 *        used a seed for the rng.
 * 
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \endclist
 *
 * Time complexity: The algorithm is a randomized las vegas algorithm. The expected
 *                  running time is O(km) where k is the value mentiened above.
 *
 */
int igraph_spanner (const igraph_t *graph,
        igraph_t *spanner,
        igraph_real_t stretch,
        igraph_vector_t *weights,
        igraph_vector_t *spanner_weight) { 

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph), edge;
    long int i = 0, nlen, neighbor, v_center;
    double sample_prob, size_limit, k = (stretch + 1) / 2, weight;
    igraph_t residual_graph;
    igraph_vector_t clustering, centers, sampled_centers, lightest_neighbor, lightest_weight, residual_weight;
    igraph_vector_t new_clustering, eids, spanner_weight_temp;
    igraph_vector_int_t *neis;
    igraph_real_t generated_number;
    igraph_integer_t from, to;
    igraph_inclist_t inclist;
    igraph_es_t es;
    igraph_heap_t edges_to_remove, edges_to_add;

    if (spanner == 0) {
        return IGRAPH_SUCCESS;
    }


    // If the weights is empty (the graph is unweigthed) then assign a weight vector of ones
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
    // Initialized the spanner and the weight vector
    IGRAPH_CHECK(igraph_empty(spanner, no_of_nodes, IGRAPH_UNDIRECTED));
    IGRAPH_VECTOR_INIT_FINALLY(&spanner_weight_temp, 0);

    // Phase 1: forming the clusters
    // Create a residual graph with V' = V and E' = E. If the graph is unweighted than each edge has
    // weight of 1.
    IGRAPH_CHECK(igraph_copy(&residual_graph, graph));
    IGRAPH_FINALLY(igraph_destroy, &residual_graph);

    // A vector which maps the nodes to clusters centers. At the beginning 
    // each node is its own cluster center.
    IGRAPH_CHECK(igraph_vector_init_seq(&clustering, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &clustering);
    // A vector with the list of the clustering centers. In the beginning of the 
    // algorithm all the nodes are centers to their clusters.
    IGRAPH_CHECK(igraph_vector_init_seq(&centers, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &centers);
    // A mapping vector which indicated the neighboring edge with the minimal weight
    // to each cluster central
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_neighbor, no_of_nodes);
    // A mapping vector which indicated the minimum weight to each neighboring cluster
    IGRAPH_VECTOR_INIT_FINALLY(&lightest_weight, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&eids, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&new_clustering, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&sampled_centers, no_of_nodes);
    IGRAPH_CHECK(igraph_heap_init(&edges_to_add, 0));
    IGRAPH_FINALLY(igraph_heap_destroy, &edges_to_add);
    IGRAPH_CHECK(igraph_heap_reserve(&edges_to_add, no_of_edges));
    IGRAPH_CHECK(igraph_heap_init(&edges_to_remove, 0));
    IGRAPH_FINALLY(igraph_heap_destroy, &edges_to_remove);
    IGRAPH_CHECK(igraph_heap_reserve(&edges_to_remove, no_of_edges));
    sample_prob = pow(no_of_nodes, -1 / k);
    size_limit = pow(no_of_nodes, 1+1/k);


    while (i < k-1) {
        //step 1: sample centers
        RNG_BEGIN();
        IGRAPH_CHECK(igraph_inclist_init(&residual_graph, &inclist, IGRAPH_OUT, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        igraph_vector_fill(&new_clustering, -1);
        igraph_vector_fill(&sampled_centers, 0);
        int size_centers = igraph_vector_size(&centers);
        for (int j = 0; j < size_centers; j++) {
            generated_number = RNG_UNIF01();
            long int center = VECTOR(centers)[j];
            if (generated_number < sample_prob) {
                VECTOR(sampled_centers)[center] = 1;
            }
        }
        RNG_END();

        //step 2 and 3
        for (long int v = 0; v < no_of_nodes; v++) {
            // If the cluster of v is sampled than continue
            v_center = VECTOR(clustering)[v];
            if (VECTOR(sampled_centers)[v_center]) {
                VECTOR(new_clustering)[v] = v_center;
                continue;
            }
            
            // If the node isn't inside any cluster than continue
            neis = igraph_inclist_get(&inclist, v);
            nlen = igraph_vector_int_size(neis);
            if (VECTOR(clustering)[v] == -1) {
                continue;
            }
            igraph_vector_fill(&lightest_neighbor, -1);
            igraph_vector_fill(&lightest_weight, INFINITY);
            
            // step 2: find the neighboring clusters and lightest edges to them
            long int lightest_sample_clustering_index = -1;
            IGRAPH_CHECK(igraph_i_lightest_edge_cluster(&residual_graph, 
                                                        &inclist, 
                                                        v, 
                                                        &lightest_sample_clustering_index,
                                                        &residual_weight,
                                                        &clustering,
                                                        &sampled_centers,
                                                        &lightest_neighbor,
                                                        &lightest_weight));

            //step 3: add edges to spanner
            if (lightest_sample_clustering_index == -1) {
                // add lightest edge which connect to each neighboring center 
                for (long int j = 0; j < no_of_nodes; j++) {
                    if (VECTOR(lightest_neighbor)[j] != -1) {
                        edge = VECTOR(lightest_neighbor)[j];
                        IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                    }
                }
                // remove all incident edges
                for (int j = 0; j < nlen; j++) {
                    edge = VECTOR(*neis)[j];
                    IGRAPH_CHECK(igraph_heap_push(&edges_to_remove, (igraph_real_t) edge));
                }
            } else {
                // add the edge connecting to the lightest sampled cluster
                edge =  VECTOR(lightest_neighbor)[lightest_sample_clustering_index];
                IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                double lightest_sample_weight = VECTOR(lightest_weight)[lightest_sample_clustering_index];
                VECTOR(new_clustering)[v] = lightest_sample_clustering_index;

                // add to the spanner light edges with weight less then the lightest_sample_weight
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
                    if ((neighbor_cluster == lightest_sample_clustering_index) || (weight < lightest_sample_weight)) {
                        edge = VECTOR(*neis)[j];
                        IGRAPH_CHECK(igraph_heap_push(&edges_to_remove, (igraph_real_t) edge));
                    }
                }
            }
        }

        // check if this itteration added to many edges to the spanner
        if (igraph_heap_size(&edges_to_add) > size_limit) {
            igraph_inclist_destroy(&inclist);
            IGRAPH_FINALLY_CLEAN(1);
            // free the heaps
            while (!igraph_heap_empty(&edges_to_add)) {
                igraph_heap_delete_top(&edges_to_add);
            }
            while (!igraph_heap_empty(&edges_to_remove)) {
                igraph_heap_delete_top(&edges_to_remove);
            }
            continue;
        }
        i = i+1;

        // Add the edges to the spanner
        long int edge_old = -1;
        while (!igraph_heap_empty(&edges_to_add)) {
            edge = (long int) igraph_heap_delete_top(&edges_to_add);
            if (edge != edge_old) {
                // only if the edge is new, meaning it wasn't inserted by the j-1 element in the list
                // then insert the edge to the spanner. 
                igraph_edge(&residual_graph, edge, &from, &to);
                IGRAPH_CHECK(igraph_add_edge(spanner, from, to));
                weight = VECTOR(residual_weight)[edge];
                IGRAPH_CHECK(igraph_vector_push_back(&spanner_weight_temp, weight));
            }
            edge_old = edge;
        }

        //delte edges from residual graph
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

        //free
        igraph_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    }

    //phase 2: vertex_clustering joining
    IGRAPH_CHECK(igraph_inclist_init(&residual_graph, &inclist, IGRAPH_OUT, IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    for (int v = 0; v < no_of_nodes; v++) {
        if (VECTOR(clustering)[v] != -1) {
            igraph_vector_fill(&lightest_neighbor, -1);
            igraph_vector_fill(&lightest_weight, INFINITY);
            IGRAPH_CHECK(igraph_i_lightest_edge_cluster(&residual_graph, 
                                                        &inclist, 
                                                        v, 
                                                        NULL,
                                                        &residual_weight,
                                                        &clustering,
                                                        NULL,
                                                        &lightest_neighbor,
                                                        &lightest_weight));
            for (int j = 0; j < no_of_nodes; j++) {
                if (VECTOR(lightest_neighbor)[j] != -1) {
                    edge = VECTOR(lightest_neighbor)[j];
                    IGRAPH_CHECK(igraph_heap_push(&edges_to_add, (igraph_real_t) edge));
                }
            }
        }
    }
    long int edge_old = -1;
    while (!igraph_heap_empty(&edges_to_add)) {
        edge = (long int)igraph_heap_delete_top(&edges_to_add);
        if (edge != edge_old) {
            // only if the edge is new, meaning it wasn't inserted by the j-1 element in the list
            // then insert the edge to the spanner. 
            igraph_edge(&residual_graph, edge, &from, &to);
            IGRAPH_CHECK(igraph_add_edge(spanner, from, to));
            weight = VECTOR(residual_weight)[edge];
            IGRAPH_CHECK(igraph_vector_push_back(&spanner_weight_temp, weight));
        }
        edge_old = edge;
    }
    if (spanner_weight) {
        IGRAPH_CHECK(igraph_vector_copy(spanner_weight, &spanner_weight_temp));
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
    igraph_vector_destroy(&sampled_centers);
    igraph_vector_destroy(&spanner_weight_temp);
    igraph_heap_destroy(&edges_to_remove);
    igraph_heap_destroy(&edges_to_add);
    igraph_vector_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(13);
    return IGRAPH_SUCCESS;
}
