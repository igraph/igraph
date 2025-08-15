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
    igraph_integer_t num_points = igraph_matrix_nrow(points);
    igraph_integer_t num_dims   = igraph_matrix_ncol(points);
    if (beta >= 1 && num_points > num_dims) { // large beta, subset of delaunay
        IGRAPH_CHECK(igraph_i_delaunay_edges(edges, points));
    } else { // small beta, not subset of delaunay, give complete graph.
        // Or delaunay not calculable due to point count
        igraph_integer_t numpoints = igraph_matrix_nrow(points);
        for (igraph_integer_t a = 0; a < numpoints - 1; a++) {
            for (igraph_integer_t b = a + 1; b < numpoints; b++) {
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, a));
                IGRAPH_CHECK(igraph_vector_int_push_back(edges, b));
            }
        }
    }
    return IGRAPH_SUCCESS;
}

igraph_real_t standard_r(igraph_real_t beta) {
    return beta * 0.5;
}

igraph_real_t small_r(igraph_real_t beta) {
    return 0.5 / beta;
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
// Test whether the two result sets have any points in common
// They each exclude their endpoints internally, so it doesn't need to be checked explicitly.
// Since a_results can't have a, it doesn't matter if b has it.
igraph_error_t intersection_predicate(igraph_bool_t *result, GraphBuildingResultSet *a_results, GraphBuildingResultSet *b_results) {
    *result = !is_overlap(a_results->neighbors, a_results->size(), b_results->neighbors, b_results->size());
    return IGRAPH_SUCCESS;
}

// Test whether the two result sets have found anything other than the endpoints being checked.
// Needs to check explicitly for the endpoints, since a_results having b would be a false positive.
igraph_error_t union_predicate(igraph_bool_t *result, GraphBuildingResultSet * a_results, GraphBuildingResultSet * b_results) {
    *result = true;
    for (igraph_integer_t i = 0; i < a_results->size(); i++) {
        if (a_results->neighbors[i] != b_results->current_vertex) {
            *result = false;
            return IGRAPH_SUCCESS;
        }
    }
    for (igraph_integer_t i = 0; i < b_results->size(); i++) {
        if (b_results->neighbors[i] != a_results->current_vertex) {
            *result = false;
            return IGRAPH_SUCCESS;
        }
    }
    return IGRAPH_SUCCESS;
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

igraph_error_t construct_lune_centres(igraph_vector_t *a_centre, igraph_vector_t *b_centre, igraph_integer_t a, igraph_integer_t b, igraph_real_t r, const igraph_matrix_t *points) {
    igraph_vector_t a_point, b_point;
    IGRAPH_VECTOR_INIT_FINALLY(&a_point, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&b_point, 0);

    IGRAPH_CHECK(igraph_matrix_get_row(points, &a_point, a));
    IGRAPH_CHECK(igraph_matrix_get_row(points, &b_point, b));
    igraph_integer_t dims = igraph_matrix_ncol(points);

    for (igraph_integer_t i = 0; i < dims; i++) {
        VECTOR(*a_centre)[i] = VECTOR(a_point)[i] + (r - 1) * (VECTOR(a_point)[i] - VECTOR(b_point)[i]);
        VECTOR(*b_centre)[i] = VECTOR(b_point)[i] + (r - 1) * (VECTOR(b_point)[i] - VECTOR(a_point)[i]);
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
        VECTOR(perp)[i] = (VECTOR(a_point)[i] - VECTOR(b_point)[i]) * sqrt(r * r - 0.25);
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

template < igraph_error_t filter(igraph_bool_t *result, kdTree < -1 > &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta) >
igraph_error_t filter_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_integer_t point_count = igraph_matrix_nrow(points);
    igraph_integer_t available_edges = igraph_vector_int_size(edges);
    igraph_integer_t added_edges = 0;
    ig_point_adaptor adaptor(points);
    igraph_integer_t dim = igraph_matrix_ncol(points);
    kdTree < -1 > tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();
    igraph_bool_t result;
    for (igraph_integer_t i = 0; i * 2  < available_edges; i++) {
        filter(&result, tree, VECTOR(*edges)[2 * i], VECTOR(*edges)[2 * i + 1], points, beta);
        if (result) {
            VECTOR(*edges)[added_edges * 2]     = VECTOR(*edges)[i * 2];
            VECTOR(*edges)[added_edges * 2 + 1] = VECTOR(*edges)[i * 2 + 1];
            added_edges += 1;
        }
    }
    IGRAPH_CHECK(igraph_vector_int_resize(edges, added_edges * 2));
    return IGRAPH_SUCCESS;
}

template <igraph_real_t r_calculation(igraph_real_t beta),
          igraph_error_t centre_positions(igraph_vector_t *a_centre, igraph_vector_t *b_centre, igraph_integer_t a, igraph_integer_t b, igraph_real_t beta, const igraph_matrix_t * points),
          igraph_error_t result_predicate(igraph_bool_t *result, GraphBuildingResultSet *a_results, GraphBuildingResultSet *b_results)>
igraph_error_t edge_is_present(igraph_bool_t *result, kdTree < -1 > &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta) {
    // position centres correctly
    igraph_integer_t dims = igraph_matrix_ncol(points);
    igraph_real_t r = r_calculation(beta);
    igraph_vector_t a_centre, b_centre;
    IGRAPH_VECTOR_INIT_FINALLY (&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY (&b_centre, dims);

    IGRAPH_CHECK(centre_positions(&a_centre, &b_centre, a, b, r,  points));

    igraph_real_t distance = get_sqr_distance(a, b, points) * r * r ; // nanoflann uses squared distances
    // beta term is used to scale to the actual circles, not just distance between points.
    GraphBuildingResultSet a_results(IGRAPH_INTEGER_MAX, distance); // TODO: use correct neighbor count
    GraphBuildingResultSet b_results(IGRAPH_INTEGER_MAX, distance);

    a_results.reset(a);
    b_results.reset(b);

    tree.findNeighbors(a_results, VECTOR(a_centre));
    tree.findNeighbors(b_results, VECTOR(b_centre));

    std::sort(a_results.neighbors.begin(), a_results.neighbors.begin() + a_results.size());
    std::sort(b_results.neighbors.begin(), b_results.neighbors.begin() + b_results.size());

    igraph_vector_destroy(&a_centre);
    igraph_vector_destroy(&b_centre);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(result_predicate(result, &a_results, &b_results));
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_lune_beta_skeleton
 * \brief computes the lune based beta skeleton of a spatial point set.
 *
 * \experimental
 *
 * This function constructs the graph corresponding to the lune based Delaunay triangulation
 * of an n-dimensional spatial point set.
 *
 * Two points are connected if a region between them whose shape is parameterized by beta
 * is free of other points.
 *
 * A larger beta results in a larger region, and a sparser graph.
 *
 * Values of beta \lt 1 are only supported in 2d, and are considerably slower.
 *
 * The gabriel graph is a special case of beta skeleton where beta = 1.
 *
 * The Relative Neighborhood graph is a special case of beta skeleton where beta = 2.
 *
 * \param graph A pointer to the graph that will be created.
 * \param points A matrix containing the points that will be used to create the graph.
 *     Each row is a point, dimensionality is inferred from the column count.
 *
 * \return Error code.
 */
igraph_error_t igraph_lune_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    if (beta >= 1) {
        IGRAPH_CHECK((filter_edges<edge_is_present<standard_r, construct_lune_centres, intersection_predicate >> (&potential_edges, points, beta)));
    } else {
        if (igraph_matrix_ncol(points) != 2) {
            IGRAPH_ERROR("Beta skeletons with beta < 1 are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
        }

        IGRAPH_CHECK((filter_edges<edge_is_present<small_r, construct_perp_centres, intersection_predicate >> (&potential_edges, points, beta)));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_circle_beta_skeleton
 * \brief Computes the circle based beta skeleton of a 2d spatial point set.
 *
 * \experimental
 *
 * This function constructs the graph corresponding to the circle based Delaunay triangulation
 * of a 2-dimensional spatial point set.
 *
 * Two points are connected if a region between them whose shape is parameterized by beta
 * is free of other points.
 *
 * A larger beta results in a larger region, and a sparser graph.
 *
 * Values of beta \lt 1 are considerably slower
 *
 * \param graph A pointer to the graph that will be created
 * \param points A Matrix containing the points that will be used to create the graph.
 *     Each row is a point.
 * \param beta A positive real value used to parameterize the graph.
 * \return Error code.
 */
igraph_error_t igraph_circle_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    if (igraph_matrix_ncol(points) != 2) {
        IGRAPH_ERROR("Circle based beta skeletons are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
    }

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    if (beta >= 1) {

        IGRAPH_CHECK((filter_edges<edge_is_present<standard_r, construct_perp_centres, union_predicate >> (&potential_edges, points, beta)));
    } else {
        IGRAPH_CHECK((filter_edges<edge_is_present<small_r, construct_perp_centres, intersection_predicate >> (&potential_edges, points, beta)));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

// derived from code by szabolcs

class BetaFinder {
    const igraph_real_t max_beta;
    const igraph_real_t tol;
    const igraph_integer_t ai, bi;
    const igraph_matrix_t *ps;
    const igraph_real_t ab2;

    double smallest_beta;
    double max_radius;

public:
    using DistanceType = igraph_real_t;
    using IndexType = igraph_integer_t;
    BetaFinder(double max_beta, double tol, igraph_integer_t v1, igraph_integer_t v2, const igraph_matrix_t *ps) :
        max_beta(max_beta), tol(tol), ai(v1), bi(v2), ps(ps),
        ab2(get_sqr_distance(ai, bi, ps)) {
        init();
    }

    void init() {
        clear();
    }

    void clear() {
        smallest_beta = INFINITY; // TODO: find proper igraph infinity
        max_radius = luneHalfHeight2(max_beta);
    }

    bool full() const {
        return true;
    }

    size_t size() const {
        return 1;
    }

    double luneHalfHeight2(double beta) const {
        if (beta == 0) {
            return 0;
        }
        return (ab2 / 4) * (2 * beta - 1);
    }

    igraph_real_t pointBeta(igraph_integer_t index) const {
        igraph_real_t ap2 = get_sqr_distance(ai, index, ps);
        igraph_real_t bp2 = get_sqr_distance(bi, index, ps);

        if (ap2 > bp2) {
            std::swap(ap2, bp2);
        }

        double denom = ab2 + ap2 - bp2;

        if (denom <= 0) {
            return std::numeric_limits<double>::infinity();
        }

        igraph_real_t beta = 2 * ap2 / denom;

        return beta < 1 + tol ? 0 : beta;
    }

    bool addPoint(igraph_real_t dist, igraph_integer_t index) {

        //mma::mout << "considering " << index << ", dist = " << dist << ", max_radius = " << max_radius << std::endl;
        if (dist < max_radius) {
            double beta = pointBeta(index);
            //mma::mout << "beta = " << beta << std::endl;

            if (beta < smallest_beta && beta < max_beta) {
                smallest_beta = beta;
                max_radius = luneHalfHeight2(beta);
            }
        }

        return true;
    }

    igraph_real_t worstDist() const {
        return max_radius;
    }
    void sort() {}
    igraph_real_t const thresholdBeta() const {
        return smallest_beta;
    }
};

igraph_error_t igraph_beta_weighted_gabriel_graph(igraph_t *graph, igraph_vector_t *edge_weights, const igraph_matrix_t *points, igraph_real_t max_beta) {
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    igraph_integer_t dim = igraph_matrix_ncol(points);
    igraph_integer_t point_count = igraph_matrix_nrow(points);
    ig_point_adaptor adaptor(points);

    IGRAPH_CHECK(igraph_i_delaunay_edges(&edges, points));
    igraph_integer_t edge_count = igraph_vector_int_size(&edges) / 2;
    IGRAPH_CHECK(igraph_vector_resize(edge_weights, edge_count));
    kdTree < -1 > tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

    igraph_vector_t midpoint;
    IGRAPH_VECTOR_INIT_FINALLY(&midpoint, dim);

    for (igraph_integer_t i = 0; i < edge_count; i++) {
        BetaFinder finder(max_beta, /*tolerance=*/0, VECTOR(edges)[2 * i], VECTOR(edges)[2 * i + 1], points);

        for (igraph_integer_t axis = 0; axis < dim; axis++) {
            VECTOR(midpoint)[axis] = 0.5 * (
                                         MATRIX(*points, VECTOR(edges)[2 * i],     axis)
                                         + MATRIX(*points, VECTOR(edges)[2 * i + 1], axis)
                                     );
        }

        tree.findNeighbors(finder, VECTOR(midpoint));
        VECTOR(*edge_weights)[i] = finder.thresholdBeta();

    }

    igraph_integer_t added_edges = 0;
    for (igraph_integer_t i = 0; i < edge_count; i++) {
        if (VECTOR(*edge_weights)[i] != 0) {
            VECTOR(*edge_weights)[added_edges] = VECTOR(*edge_weights)[i];
            VECTOR(edges)[added_edges * 2] = VECTOR(edges)[i * 2];
            VECTOR(edges)[added_edges * 2 + 1] = VECTOR(edges)[i * 2 + 1];
            added_edges += 1;
        }
    }

    IGRAPH_CHECK(igraph_vector_int_resize(&edges, added_edges * 2));
    IGRAPH_CHECK(igraph_vector_resize(edge_weights, added_edges));

    IGRAPH_CHECK(igraph_create(graph, &edges, point_count, false));

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&midpoint);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
