/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_matrix.h"
#include "igraph_interface.h"
#include "igraph_umap.h"
#include "igraph_constructors.h"

static igraph_error_t igraph_get_k_smallest(igraph_vector_t *values,
        igraph_vector_int_t *smallest, igraph_integer_t k, igraph_integer_t starting_index)
{
    return (IGRAPH_SUCCESS);
}

/*This creates a k-nearest neighbors graph with their distances.*/
/*first trying everything with simple slow algorithms */
static igraph_error_t igraph_knng_from_euclid(igraph_matrix_t *data,
        igraph_t *knn_graph, igraph_vector_t *knn_distances) {
    igraph_integer_t k = 10;
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(data);
    igraph_vector_t node_distances;
    igraph_vector_t data_row_i;
    igraph_vector_t data_row_j;
    igraph_vector_int_t closest;
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INIT_FINALLY(&node_distances, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&data_row_i, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&data_row_j, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&closest, k);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_matrix_get_row(data, &data_row_i, i);
        for (igraph_integer_t j = i; j < no_of_nodes; j++) {
            igraph_matrix_get_row(data, &data_row_j, j);
            igraph_vector_sub(&data_row_j, &data_row_i);
            VECTOR(node_distances)[j] = igraph_vector_sumsq(&data_row_j);
        }
        igraph_get_k_smallest(&node_distances, &closest, k, i); //uses i as a starting point. Very ugly interface
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(&closest); j++) { //size is usually k, but not for last few nodes
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, VECTOR(closest)[j]));
            IGRAPH_CHECK(igraph_vector_push_back(knn_distances, VECTOR(node_distances)[VECTOR(closest)[j]]));
        }


    }
    IGRAPH_CHECK(igraph_create(knn_graph, &edges, no_of_nodes, IGRAPH_UNDIRECTED));
    igraph_vector_destroy(&node_distances);
    igraph_vector_destroy(&data_row_i);
    igraph_vector_destroy(&data_row_j);
    igraph_vector_int_destroy(&closest);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(5);
    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_umap_find_open_sets(igraph_t *knn_graph,
        igraph_vector_t *distances, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {

    return (IGRAPH_SUCCESS);
}


static igraph_error_t igraph_umap_edge_weights(igraph_matrix_t *data,
        igraph_t *umap_graph, igraph_vector_t *open_set_size,
        igraph_vector_t *open_set_decay) {

    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_umap_layout(igraph_t *umap_graph,
        igraph_matrix_t *layout) {

    return (IGRAPH_SUCCESS);
}
igraph_error_t igraph_umap(igraph_matrix_t *data, igraph_matrix_t *layout) {
    igraph_t knn_graph;
    igraph_t umap_graph;
    igraph_vector_t open_set_sizes;
    igraph_vector_t open_set_decays;
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(data);
    igraph_vector_t knn_distances;

    IGRAPH_VECTOR_INIT_FINALLY(&open_set_sizes, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&open_set_decays, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&knn_distances, 0);
    IGRAPH_CHECK(igraph_matrix_resize(layout, no_of_nodes, 2));

    IGRAPH_CHECK(igraph_knng_from_euclid(data, &knn_graph, &knn_distances));
    IGRAPH_CHECK(igraph_umap_find_open_sets(&knn_graph, &knn_distances, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_edge_weights(data, &umap_graph, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_layout(&umap_graph, layout));

    igraph_vector_destroy(&open_set_sizes);
    igraph_vector_destroy(&open_set_decays);
    igraph_vector_destroy(&knn_distances);
    IGRAPH_FINALLY_CLEAN(3);
    return (IGRAPH_SUCCESS);
}
