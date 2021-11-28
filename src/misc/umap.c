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

static igraph_error_t igraph_knng_from_euclid(igraph_matrix_t *data,
        igraph_t *knn_graph, igraph_vector_t *knn_distances) {


    return (IGRAPH_SUCCESS);
}

static igraph_error_t igraph_umap_find_open_sets(igraph_t *knn_graph,
        igraph_vector_t *distances, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {

    return (IGRAPH_SUCCESS);
}


igraph_error_t igraph_umap_edge_weights(igraph_matrix_t *data,
        igraph_t *umap_graph, igraph_vector_t *open_set_size,
        igraph_vector_t *open_set_decay) {

    return (IGRAPH_SUCCESS);
}

igraph_error_t igraph_umap_layout(igraph_t *umap_graph,
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
