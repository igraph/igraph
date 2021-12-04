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
#include "igraph_layout.h"

/*open set size is just the size of the distance to the closest neighbor*/
/*the decay depends on the rest of the neighbors */

/*now to calculate the decay*/
/*we could curve fit an offset exponential.*/
/*for now we can do somehting simpler which doesn't use all the k points*/

/*this picks the largest distance and calculates the decay from that. Just a basic experiment*/
static igraph_error_t igraph_umap_find_open_sets(igraph_t *knn_graph,
        igraph_vector_t *distances, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {
    igraph_integer_t no_of_nodes = igraph_vcount(knn_graph);
    igraph_vector_int_t eids;
    igraph_real_t largest;

    /*this picks the shortest distance*/
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_incident(knn_graph, &eids, i, IGRAPH_ALL);
        VECTOR(*open_set_sizes)[i] = VECTOR(*distances)[VECTOR(eids)[0]];
        largest = -1;
        for (igraph_integer_t j = 1; j < igraph_vector_int_size(&eids); j++) {
            if (VECTOR(*distances)[VECTOR(eids)[j]] < VECTOR(*open_set_sizes)[i]) {
                VECTOR(*open_set_sizes)[i] = VECTOR(*distances)[VECTOR(eids)[j]];
            }
            if (VECTOR(*distances)[VECTOR(eids)[j]] > largest) {
                largest = VECTOR(*distances)[VECTOR(eids)[j]];
            }
        }
        VECTOR(*open_set_decays)[i] = (largest - VECTOR(*open_set_sizes)[i]);

    }
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(1);
    return (IGRAPH_SUCCESS);
}


igraph_error_t igraph_umap_decay(igraph_real_t *weight, igraph_real_t distance, igraph_real_t open_set_size, igraph_real_t open_set_decay)
{
    return (IGRAPH_SUCCESS);
}


/*a point within the open_set_size gets distance of 1, they can't be actually in it....they're also not really 'open' disks btw.*/
/* and then 1 + (distance - open_set_size) / (e ^ -(distance - open_set_size) * decay)? */
/* which is 1 + (distance - open_set_size) * (e ^ (distance - open_set_size) * decay)? */
/* just start with a linear decay */
static igraph_error_t igraph_umap_edge_weights(igraph_matrix_t *data,
        igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {

    /* we go over all the nodes, and then over all the nodes, and then add an edge with the weight dependent on the open set and decay*/
    /* and no edge if this weight would be 0 or lower */
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(data);
    igraph_vector_int_t edges;
    igraph_vector_t data_row_i;
    igraph_vector_t data_row_j;
    igraph_real_t distance;
    igraph_real_t weight_a;
    igraph_real_t weight_b;
    igraph_real_t weight;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_matrix_get_row(data, &data_row_i, i);
        for (igraph_integer_t j = i; j < no_of_nodes; j++) {
            igraph_matrix_get_row(data, &data_row_j, j);
            igraph_vector_sub(&data_row_j, &data_row_i);
            distance = igraph_vector_sumsq(&data_row_j);
            IGRAPH_CHECK(igraph_umap_decay(&weight_a,  distance, VECTOR(*open_set_sizes)[i], VECTOR(*open_set_decays)[i]));
            IGRAPH_CHECK(igraph_umap_decay(&weight_b, distance, VECTOR(*open_set_sizes)[j], VECTOR(*open_set_decays)[j]));
            if (weight_a <= 0 || weight_b <= 0)
                continue;
            weight = weight_a + weight_b - weight_a * weight_b;
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
            IGRAPH_CHECK(igraph_vector_push_back(umap_weights, weight));
        }
    }

    IGRAPH_CHECK(igraph_create(umap_graph, &edges, no_of_nodes, IGRAPH_UNDIRECTED));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_get_gradient(igraph_matrix_t *gradient, igraph_matrix_t *layout, igraph_t *umap_graph, igraph_vector_t *umap_weights)
{
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_vector_int_t eids;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_incident(umap_graph, &eids, i, IGRAPH_ALL);
        MATRIX(*gradient, i, 0) = 0;
        MATRIX(*gradient, i, 1) = 0;
        for (igraph_integer_t j = 1; j < igraph_vector_int_size(&eids); j++) {
            igraph_integer_t eid = VECTOR(eids)[j];
            igraph_real_t x = MATRIX(*layout, i, 0);
            igraph_real_t y = MATRIX(*layout, i, 1);
            igraph_integer_t other = IGRAPH_OTHER(umap_graph, eid, i);
            igraph_real_t other_x = MATRIX(*layout, other, 0) ;
            igraph_real_t other_y = MATRIX(*layout, other, 1) ;
            igraph_real_t x_diff = (x - other_x);
            igraph_real_t y_diff = (y - other_y);
            igraph_real_t distance = (x_diff * x_diff + y_diff * y_diff);
            igraph_real_t weight_diff = -VECTOR(*umap_weights)[j] - 1 / distance;
            igraph_real_t gradient_x = weight_diff * x_diff / distance;
            igraph_real_t gradient_y = weight_diff * y_diff / distance;
            MATRIX(*gradient, i, 0) += gradient_x;
            MATRIX(*gradient, i, 1) += gradient_y;
        }
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_umap_layout(igraph_t *umap_graph, igraph_vector_t *umap_weights,
        igraph_matrix_t *layout) {
    igraph_integer_t epochs = 100;
    igraph_real_t learning_rate = 0.1;
    igraph_matrix_t gradient;
    igraph_matrix_init(&gradient, igraph_matrix_nrow(layout), igraph_matrix_ncol(layout));

    igraph_layout_random(umap_graph, layout);
    for (igraph_integer_t e = 0; e < epochs; e++) {
        igraph_get_gradient(&gradient, layout, umap_graph, umap_weights);
        igraph_matrix_scale(&gradient, learning_rate);
        igraph_matrix_sub(layout, &gradient);
    }

    IGRAPH_FINALLY_CLEAN(1);

    return (IGRAPH_SUCCESS);
}

igraph_error_t igraph_layout_umap(igraph_t *graph, igraph_vector_t *weights, igraph_matrix_t *layout) {
    igraph_vector_t open_set_sizes;
    igraph_vector_t open_set_decays;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_t knn_distances;
    igraph_vector_t umap_weights;

    IGRAPH_VECTOR_INIT_FINALLY(&open_set_sizes, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&open_set_decays, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&knn_distances, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&umap_weights, 0);
    IGRAPH_CHECK(igraph_matrix_resize(layout, no_of_nodes, 2));

    IGRAPH_CHECK(igraph_umap_find_open_sets(graph, weights, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_edge_weights(data, &umap_graph, &umap_weights, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_layout(&umap_graph, &umap_weights, layout));

    igraph_vector_destroy(&open_set_sizes);
    igraph_vector_destroy(&open_set_decays);
    igraph_vector_destroy(&knn_distances);
    igraph_vector_destroy(&umap_weights);
    IGRAPH_FINALLY_CLEAN(4);
    return (IGRAPH_SUCCESS);
}
