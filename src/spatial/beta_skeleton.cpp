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

#include "igraph_constructors.h"
#include "igraph_matrix.h"
#include "igraph_spatial.h"

#include "igraph_vector.h"

#include "spatial/spatial_internal.h"
#include "spatial/nanoflann_internal.hpp"
#include "nanoflann/nanoflann.hpp"

#include "igraph_error.h"
#include <algorithm>
#include <cstdlib>
#include <vector>

template <igraph_integer_t Dimension>
using kdTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<igraph_real_t, ig_point_adaptor>, ig_point_adaptor, Dimension, igraph_integer_t>;

kdTree<0> build_tree() {

}

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

igraph_error_t filter_small_edges(igraph_vector_int_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    return IGRAPH_SUCCESS;
}

igraph_bool_t is_overlap(std::vector<igraph_integer_t> &a, std::vector<igraph_integer_t> &b) {
    igraph_integer_t
        a_pos = 0,
        a_size = a.size(),
        b_pos = 0,
        b_size = b.size();

    while (a_pos < a_size && b_pos < b_size) {
        if (a[a_pos] == b[b_pos]) return true;
        if (a[a_pos] < b[b_pos]) {
            a_pos++;
        } else {
            b_pos++;
        }
    }
    return false;
}

igraph_real_t sqrDistance(igraph_vector_t *a, igraph_vector_t *b) {
    igraph_integer_t size = igraph_vector_size(a);
    igraph_real_t accumulator = 0;
    igraph_real_t temp;
    for (igraph_integer_t i = 0; i < size; i++) {
        temp = abs(VECTOR(*a)[i] - VECTOR(*b)[i]);
        accumulator += temp * temp;
    }
    return accumulator;
}

igraph_bool_t lune_is_present(kdTree<0> &tree, igraph_integer_t a, igraph_integer_t b, igraph_matrix_t *points, igraph_real_t beta) {

    igraph_vector_t a_point, b_point;

    IGRAPH_CHECK(igraph_matrix_get_row(points, &a_point, a));
    IGRAPH_CHECK(igraph_matrix_get_row(points, &b_point, b));

    igraph_vector_t a_centre, b_centre;


    igraph_integer_t dims = igraph_matrix_ncol(points);

    IGRAPH_VECTOR_INIT_FINALLY(&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY(&b_centre, dims);

    for (igraph_integer_t i = 0; i < dims; i++) {
        VECTOR(a_centre)[i] = VECTOR(a_point)[i] + (0.5*beta - 1) * (VECTOR(a_centre)[i] - VECTOR(b_centre)[i]);
        VECTOR(b_centre)[i] = VECTOR(b_point)[i] + (0.5*beta - 1) * (VECTOR(b_centre)[i] - VECTOR(a_centre)[i]);
    }

    igraph_integer_t distance = sqrDistance(&a_point, &b_point); // nanoflann uses squared distances

    GraphBuildingResultSet a_results(10000, distance);
    GraphBuildingResultSet b_results(10000, distance);

    a_results.reset(-1);
    b_results.reset(-1);

    tree.findNeighbors(a_results, VECTOR(a_centre));
    tree.findNeighbors(b_results, VECTOR(b_centre));

    std::sort(a_results.neighbors.begin(), a_results.neighbors.end());
    std::sort(b_results.neighbors.begin(), b_results.neighbors.end());

    return is_overlap(a_results.neighbors, b_results.neighbors);
}

igraph_error_t filter_lune_edges(igraph_vector_int_t *edges, igraph_matrix_t *points, igraph_real_t beta) {
    if (beta < 1) {
        IGRAPH_CHECK(filter_small_edges(edges, points, beta));
        return IGRAPH_SUCCESS;
    }

    igraph_integer_t point_count = igraph_matrix_nrow(points);
    igraph_integer_t added_edges = 0;
    ig_point_adaptor adaptor(points);
    igraph_integer_t dim = igraph_matrix_ncol(points);
    kdTree<0> tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));

    for (igraph_integer_t i = 0; i < point_count; i++) {
        if (lune_is_present(tree, VECTOR(*edges)[2*i], VECTOR(*edges)[2*i+1], points, beta)) {
            VECTOR(*edges)[added_edges * 2]     = VECTOR(*edges)[i*2];
            VECTOR(*edges)[added_edges * 2 + 1] = VECTOR(*edges)[i*2 + 1];
        }
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_lune_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));

    IGRAPH_CHECK(filter_lune_edges(&potential_edges, points, beta));

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, false, false));

    return IGRAPH_SUCCESS;
}


igraph_error_t igraph_circle_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {

    return IGRAPH_SUCCESS;
}

