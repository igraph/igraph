/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_constants.h"
#include "igraph_constructors.h"
#include "igraph_matrix.h"
#include "igraph_spatial.h"

#include "igraph_types.h"
#include "igraph_vector.h"

#include "spatial/spatial_internal.h"
#include "spatial/nanoflann_internal.hpp"
#include "nanoflann/nanoflann.hpp"

#include "igraph_error.h"
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <vector>



template <igraph_integer_t Dimension>
using kdTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<igraph_real_t, ig_point_adaptor>, ig_point_adaptor, Dimension, igraph_integer_t>;

igraph_error_t beta_skeleton_edge_superset(igraph_vector_int_t *edges, const igraph_matrix_t *points, igraph_real_t beta) {
    if (beta >= 1) { // large beta, subset of delaunay
        IGRAPH_CHECK(igraph_i_delaunay_edges(edges, points));
    } else { // small beta, not subset of delaunay, give complete graph.
        igraph_integer_t numpoints = igraph_matrix_nrow(points);
        for (igraph_integer_t a = 0; a < numpoints -1; a++) {
            for (igraph_integer_t b = a+1; b < numpoints; b++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, a));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, b));
            }
        }
    }
    return IGRAPH_SUCCESS;
}


igraph_bool_t is_overlap(std::vector<igraph_integer_t> &a, igraph_integer_t a_size, std::vector<igraph_integer_t> &b, igraph_integer_t b_size) {
    igraph_integer_t
        a_pos = 0,
        b_pos = 0;
    while (a_pos < a_size && b_pos < b_size) {
        if (a[a_pos] == b[b_pos]) {
            return true;
        }
        if (a[a_pos] < b[b_pos]) {
            a_pos++;
        } else {
            b_pos++;
        }
    }
    return false;
}

igraph_real_t sqr_distance(igraph_vector_t *a, igraph_vector_t *b) {
    igraph_integer_t size = igraph_vector_size(a);
    igraph_real_t accumulator = 0;
    igraph_real_t temp;
    for (igraph_integer_t i = 0; i < size; i++) {
        temp = abs(VECTOR(*a)[i] - VECTOR(*b)[i]);
        accumulator += temp * temp;
    }
    return accumulator;
}

igraph_real_t get_sqr_distance(igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points) {
    igraph_real_t distance = 0;
    igraph_real_t temp;
    igraph_integer_t dims = igraph_matrix_ncol(points);
    for (igraph_integer_t i = 0; i < dims; i++) {
        temp = abs(MATRIX(*points, a, i) - MATRIX(*points, b, i));
        distance += temp * temp;
    }
    return distance;
}

igraph_error_t construct_lune_centres(igraph_vector_t *a_centre, igraph_vector_t *b_centre, igraph_integer_t a, igraph_integer_t b, igraph_real_t beta, const igraph_matrix_t *points){
    igraph_vector_t a_point, b_point;
    IGRAPH_VECTOR_INIT_FINALLY(&a_point, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&b_point, 0);

    IGRAPH_CHECK(igraph_matrix_get_row(points, &a_point, a));
    IGRAPH_CHECK(igraph_matrix_get_row(points, &b_point, b));
    igraph_integer_t dims = igraph_matrix_ncol(points);

    for (igraph_integer_t i = 0; i < dims; i++) {
        VECTOR(*a_centre)[i] = VECTOR(a_point)[i] + (0.5*beta - 1) * (VECTOR(a_point)[i] - VECTOR(b_point)[i]);
        VECTOR(*b_centre)[i] = VECTOR(b_point)[i] + (0.5*beta - 1) * (VECTOR(b_point)[i] - VECTOR(a_point)[i]);
    }

    igraph_vector_destroy(&a_point);
    igraph_vector_destroy(&b_point);

    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

igraph_error_t construct_perp_centres(igraph_vector_t *a_centre, igraph_vector_t *b_centre, igraph_integer_t a, igraph_integer_t b, igraph_real_t r, const igraph_matrix_t *points) {
 igraph_vector_t a_point, b_point;
    IGRAPH_VECTOR_INIT_FINALLY(&a_point, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&b_point, 0);

    IGRAPH_CHECK(igraph_matrix_get_row(points, &a_point, a));
    IGRAPH_CHECK(igraph_matrix_get_row(points, &b_point, b));


    igraph_integer_t dims = igraph_matrix_ncol(points);

    igraph_vector_t mid, perp;

    IGRAPH_VECTOR_INIT_FINALLY(&mid, dims);
    IGRAPH_VECTOR_INIT_FINALLY(&perp, dims);
    for (igraph_integer_t i = 0; i < dims; i++) {
        VECTOR(mid)[i] = (VECTOR(a_point)[i] + VECTOR(b_point)[i]) * 0.5;
        VECTOR(perp)[i] = (VECTOR(a_point)[i] - VECTOR(b_point)[i]) * sqrt(r*r -0.25);
    }

    // Since this is only well defined for 2d, a manual 90 degree rotation works and is simpler.
    // the rotation being x = -y, y = x, 90 degrees counter-clockwise.
    igraph_real_t temp = VECTOR(perp)[0];
    VECTOR(perp)[0] = - VECTOR(perp)[1];
    VECTOR(perp)[1] = temp;

    for (igraph_integer_t i = 0; i < dims; i++) {
        VECTOR(*a_centre)[i] = VECTOR(mid)[i] + VECTOR(perp)[i];
        VECTOR(*b_centre)[i] = VECTOR(mid)[i] - VECTOR(perp)[i];
    }

    igraph_vector_destroy(&mid);
    igraph_vector_destroy(&perp);
    igraph_vector_destroy(&a_point);
    igraph_vector_destroy(&b_point);
    IGRAPH_FINALLY_CLEAN(4);
    return IGRAPH_SUCCESS;
}

igraph_error_t small_is_present(igraph_bool_t *result, kdTree<-1> &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta) {
        // position centres correctly
    igraph_real_t r = 0.5/beta;

    igraph_integer_t dims = igraph_matrix_ncol(points);

    igraph_vector_t a_centre, b_centre;

    IGRAPH_VECTOR_INIT_FINALLY(&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY(&b_centre, dims);

    IGRAPH_CHECK(construct_perp_centres(&a_centre, &b_centre, a, b, r, points));

    igraph_real_t distance = get_sqr_distance(a, b, points) * r * r ; // nanoflann uses squared distances
    // beta term is used to scale to the actual circles, not just distance between points.
    GraphBuildingResultSet a_results(IGRAPH_INTEGER_MAX, distance); // TODO: use correct neighbor count
    GraphBuildingResultSet b_results(IGRAPH_INTEGER_MAX, distance);

    a_results.reset(a);
    b_results.reset(b);

    tree.findNeighbors(a_results, VECTOR(a_centre));
    tree.findNeighbors(b_results, VECTOR(b_centre));

    std::sort(a_results.neighbors.begin(), a_results.neighbors.begin()+a_results.size());
    std::sort(b_results.neighbors.begin(), b_results.neighbors.begin()+b_results.size());

    igraph_vector_destroy(&a_centre);
    igraph_vector_destroy(&b_centre);
    IGRAPH_FINALLY_CLEAN(2);
    *result = !is_overlap(a_results.neighbors, a_results.size(), b_results.neighbors, b_results.size());
    return IGRAPH_SUCCESS;
}

template <igraph_error_t filter(igraph_bool_t *result, kdTree<-1> &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta)>
igraph_error_t filter_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points, igraph_real_t beta) {
    if (igraph_matrix_ncol(points) != 2) {
        IGRAPH_ERROR("Beta-skeletons with beta < 1 are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
    }
     igraph_integer_t point_count = igraph_matrix_nrow(points);
    igraph_integer_t available_edges = igraph_vector_int_size(edges);
    igraph_integer_t added_edges = 0;
    ig_point_adaptor adaptor(points);
    igraph_integer_t dim = igraph_matrix_ncol(points);
    kdTree<-1> tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();
    igraph_bool_t result;
    for (igraph_integer_t i = 0; i * 2  < available_edges; i++) {
        filter(&result, tree, VECTOR(*edges)[2*i], VECTOR(*edges)[2*i+1], points, beta);
        if (result) {
            VECTOR(*edges)[added_edges * 2]     = VECTOR(*edges)[i*2];
            VECTOR(*edges)[added_edges * 2 + 1] = VECTOR(*edges)[i*2 + 1];
            added_edges+=1;
        }
    }
    IGRAPH_CHECK(igraph_vector_int_resize(edges, added_edges * 2));
    return IGRAPH_SUCCESS;
}

igraph_error_t circle_is_present(igraph_bool_t *result, kdTree<-1> &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta) {
        // position centres correctly
    igraph_vector_t a_centre, b_centre;
    igraph_real_t r = 0.5 * beta;
    igraph_integer_t dims = igraph_matrix_ncol(points);

    IGRAPH_VECTOR_INIT_FINALLY(&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY(&b_centre, dims);

    IGRAPH_CHECK(construct_perp_centres(&a_centre, &b_centre, a, b, r, points));

    igraph_real_t distance = get_sqr_distance(a, b, points) * r * r ; // nanoflann uses squared distances
    // beta term is used to scale to the actual circles, not just distance between points.
    GraphBuildingResultSet a_results(IGRAPH_INTEGER_MAX, distance); // TODO: use correct neighbor count
    GraphBuildingResultSet b_results(IGRAPH_INTEGER_MAX, distance);

    a_results.reset(a);
    b_results.reset(b);

    tree.findNeighbors(a_results, VECTOR(a_centre));
    tree.findNeighbors(b_results, VECTOR(b_centre));

    std::sort(a_results.neighbors.begin(), a_results.neighbors.begin()+a_results.size());
    std::sort(b_results.neighbors.begin(), b_results.neighbors.begin()+b_results.size());

    igraph_vector_destroy(&a_centre);
    igraph_vector_destroy(&b_centre);
    IGRAPH_FINALLY_CLEAN(2);

    *result = true;
    for (igraph_integer_t i = 0; i < a_results.size(); i++) {
        if (a_results.neighbors[i] != b){
            *result = false;
            return IGRAPH_SUCCESS;
        }
    }
    for (igraph_integer_t i = 0; i < b_results.size(); i++) {
        if (b_results.neighbors[i] != a) {
            *result = false;
            return IGRAPH_SUCCESS;
        }
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t lune_is_present(igraph_bool_t *result, kdTree<-1> &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta) {
    // position centres correctly

    igraph_integer_t dims = igraph_matrix_ncol(points);

    igraph_vector_t a_centre, b_centre;
    IGRAPH_VECTOR_INIT_FINALLY (&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY (&b_centre, dims);

    IGRAPH_CHECK(construct_lune_centres(&a_centre, &b_centre, a, b, beta,  points));

    igraph_real_t distance = get_sqr_distance(a,b, points) * (0.5*beta) * (0.5*beta) ; // nanoflann uses squared distances
    // beta term is used to scale to the actual circles, not just distance between points.
    GraphBuildingResultSet a_results(IGRAPH_INTEGER_MAX, distance); // TODO: use correct neighbor count
    GraphBuildingResultSet b_results(IGRAPH_INTEGER_MAX, distance);

    a_results.reset(a);
    b_results.reset(b);

    tree.findNeighbors(a_results, VECTOR(a_centre));
    tree.findNeighbors(b_results, VECTOR(b_centre));

    std::sort(a_results.neighbors.begin(), a_results.neighbors.begin()+a_results.size());
    std::sort(b_results.neighbors.begin(), b_results.neighbors.begin()+b_results.size());

    igraph_vector_destroy(&a_centre);
    igraph_vector_destroy(&b_centre);
    IGRAPH_FINALLY_CLEAN(2);

    *result = !is_overlap(a_results.neighbors, a_results.size(), b_results.neighbors, b_results.size());
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_lune_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    if (beta >= 1) {
        IGRAPH_CHECK(filter_edges<lune_is_present>(&potential_edges, points, beta));
    } else {
        IGRAPH_CHECK(filter_edges<small_is_present>(&potential_edges, points, beta));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, false, false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


igraph_error_t igraph_circle_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    if (igraph_matrix_ncol(points) != 2) {
        IGRAPH_ERROR("Circle based beta skeletons are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
    }

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    if (beta >= 1) {
        IGRAPH_CHECK(filter_edges<circle_is_present>(&potential_edges, points, beta));
    } else {
        IGRAPH_CHECK(filter_edges<small_is_present>(&potential_edges, points, beta));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, false, false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

