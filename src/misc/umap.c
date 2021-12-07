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
        VECTOR(*open_set_decays)[i] = (largest - VECTOR(*open_set_sizes)[i]); /*TODO find a sensible decay, we could actually fit a function like we're supposed to */

    }
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(1);
    return (IGRAPH_SUCCESS);
}


igraph_error_t igraph_umap_decay(igraph_real_t *probability, igraph_real_t distance, igraph_real_t open_set_size, igraph_real_t open_set_decay)
{
    /* TODO: implement. this should use the fuzzy ball and distance to find
     * the edge weight between two points, which should be between 0 and 1*/
    if (distance < open_set_size) {
        *probability = .99;
        return (IGRAPH_SUCCESS);
    }
    *probability = .99 - (distance - open_set_size) * open_set_decay;
    if (distance < 0.) {
        *probability = 0.;
    }
    return (IGRAPH_SUCCESS);
}


/*a point within the open_set_size gets distance of 1, they can't be actually in it....they're also not really 'open' disks btw.*/
/* and then 1 + (distance - open_set_size) / (e ^ -(distance - open_set_size) * decay)? */
/* which is 1 + (distance - open_set_size) * (e ^ (distance - open_set_size) * decay)? */
/* just start with a linear decay */
/* TODO: don't use raw data, only knn graph */
static igraph_error_t igraph_umap_edge_weights(igraph_t *graph, igraph_vector_t *distances,
        igraph_vector_t *umap_weights, igraph_vector_t *open_set_sizes,
        igraph_vector_t *open_set_decays) {

    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t weight;
    igraph_real_t weight_previous;

    igraph_vector_resize(umap_weights, igraph_vector_size(distances));
    igraph_vector_null(umap_weights);
    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        IGRAPH_CHECK(igraph_umap_decay(&weight,  VECTOR(*distances)[i],  VECTOR(*open_set_sizes)[i], VECTOR(*open_set_decays)[i]));
        weight_previous = VECTOR(*umap_weights)[i];
        if (weight_previous > 0)
            weight = weight + weight_previous - weight * weight_previous;
        if (weight > 0)
            VECTOR(*umap_weights)[i] = weight;
        else
            VECTOR(*umap_weights)[i] = 0;
    }

    return IGRAPH_SUCCESS;
}

/*Gives the partial derivative to with respect to b */
typedef igraph_error_t igraph_partial_derivative_2d(igraph_real_t a, igraph_real_t b, igraph_real_t *derivative);

static igraph_error_t igraph_cross_entropy_derivative(igraph_real_t a, igraph_real_t b, igraph_real_t *derivative)
{
    *derivative = (b - a) / (b * (1 - b));
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_get_gradient(igraph_matrix_t *gradient, igraph_matrix_t *layout, igraph_t *umap_graph, igraph_vector_t *umap_weights, igraph_partial_derivative_2d grad)
{
    /*TODO use cross entropy, */
    igraph_integer_t no_of_nodes = igraph_matrix_nrow(layout);
    igraph_vector_int_t eids;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&eids, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_incident(umap_graph, &eids, i, IGRAPH_ALL);
        MATRIX(*gradient, i, 0) = 0;
        MATRIX(*gradient, i, 1) = 0;
        for (igraph_integer_t j = 0; j < igraph_vector_int_size(&eids); j++) {
            igraph_integer_t eid = VECTOR(eids)[j];
            igraph_real_t x = MATRIX(*layout, i, 0);
            igraph_real_t y = MATRIX(*layout, i, 1);
            igraph_integer_t other = IGRAPH_OTHER(umap_graph, eid, i);
            igraph_real_t other_x = MATRIX(*layout, other, 0) ;
            igraph_real_t other_y = MATRIX(*layout, other, 1) ;
            igraph_real_t x_diff = (x - other_x);
            igraph_real_t y_diff = (y - other_y);
            igraph_real_t distance = (x_diff * x_diff + y_diff * y_diff);
            igraph_real_t weight_2d = 1 / distance;
            igraph_real_t d;
            IGRAPH_CHECK(grad(weight_2d, VECTOR(*umap_weights)[j], &d));
            igraph_real_t gradient_x = d * x_diff / distance;
            igraph_real_t gradient_y = d * y_diff / distance;
            MATRIX(*gradient, i, 0) += gradient_x;
            MATRIX(*gradient, i, 1) += gradient_y;
        }
    }
    igraph_vector_int_destroy(&eids);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_umap_layout(igraph_t *umap_graph, igraph_vector_t *umap_weights,
        igraph_matrix_t *layout) {
    igraph_integer_t epochs = 5000;
    igraph_real_t learning_rate = 0.001;
    igraph_matrix_t gradient;
    igraph_layout_random(umap_graph, layout);
    IGRAPH_MATRIX_INIT_FINALLY(&gradient, igraph_matrix_nrow(layout), igraph_matrix_ncol(layout));

    for (igraph_integer_t e = 0; e < epochs; e++) {
        igraph_get_gradient(&gradient, layout, umap_graph, umap_weights, igraph_cross_entropy_derivative);
        //printf("gradient:\n");
        //igraph_matrix_print(&gradient);
        //printf("\n");
        igraph_matrix_scale(&gradient, learning_rate);
        igraph_matrix_sub(layout, &gradient);
        //printf("layout:\n");
        //igraph_matrix_print(layout);
        //printf("\n");
    }

    igraph_matrix_destroy(&gradient);
    IGRAPH_FINALLY_CLEAN(1);
    return (IGRAPH_SUCCESS);
}

igraph_error_t igraph_layout_umap(igraph_t *graph, igraph_vector_t *distances, igraph_matrix_t *layout) {
    igraph_vector_t open_set_sizes;
    igraph_vector_t open_set_decays;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_t umap_weights;

    IGRAPH_VECTOR_INIT_FINALLY(&open_set_sizes, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&open_set_decays, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&umap_weights, 0);

    IGRAPH_CHECK(igraph_umap_find_open_sets(graph, distances, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_edge_weights(graph, distances, &umap_weights, &open_set_sizes,
                &open_set_decays));
    IGRAPH_CHECK(igraph_umap_layout(graph, &umap_weights, layout));

    igraph_vector_destroy(&open_set_sizes);
    igraph_vector_destroy(&open_set_decays);
    igraph_vector_destroy(&umap_weights);
    IGRAPH_FINALLY_CLEAN(3);
    return (IGRAPH_SUCCESS);
}
