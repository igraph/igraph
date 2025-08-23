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
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_spatial.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "spatial/spatial_internal.h"
#include "spatial/nanoflann_internal.hpp"

#include <cfloat>

#define TOLERANCE (128 * DBL_EPSILON)

// Some methods to get distances between two vectors, including when one or both are embedded in a matrix.
static inline igraph_real_t ind_ind_sqr_distance(igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points) {
    igraph_real_t distance = 0;
    igraph_real_t temp;
    igraph_integer_t dims = igraph_matrix_ncol(points);
    for (igraph_integer_t i = 0; i < dims; i++) {
        temp = MATRIX(*points, a, i) - MATRIX(*points, b, i);
        distance += temp * temp;
    }
    return distance;
}

static inline igraph_real_t vec_vec_sqr_dist(const igraph_vector_t *a, const igraph_vector_t *b) {
    igraph_real_t distance = 0;
    igraph_real_t temp;
    igraph_integer_t dims = igraph_vector_size(a);
    for (igraph_integer_t i = 0; i < dims; i++) {
        temp = VECTOR(*a)[i] - VECTOR(*b)[i];
        distance += temp * temp;
    }
    return distance;
}

static inline igraph_real_t vec_ind_sqr_dist(const igraph_vector_t *a, const igraph_integer_t b, const igraph_matrix_t *points) {
    igraph_real_t distance = 0;
    igraph_real_t temp;
    igraph_integer_t dims = igraph_matrix_ncol(points);
    for (igraph_integer_t i = 0; i < dims; i++) {
        temp = VECTOR(*a)[i] - MATRIX(*points, b, i);
        distance += temp * temp;
    }
    return distance;
}

// Adapted from code by Szabolcs at https://github.com/szhorvat/IGraphM
// Used in is_union_empty
class NeighborCounts {
    const igraph_real_t radius; // L2 search radius
    const igraph_integer_t a, b; // edge endpoits; excluded from the count
    igraph_bool_t short_circuit;
    igraph_integer_t found;

public:
    // Boilerplate for nanoflann
    using DistanceType = igraph_real_t;
    NeighborCounts(igraph_real_t radius, igraph_integer_t a, igraph_integer_t b, igraph_bool_t short_circuit)
        : radius(radius), a(a), b(b), short_circuit(short_circuit) {
        init();
    }

    void init() {
        clear();
    }
    void clear() {
        found = 0;
    }
    size_t size() const {
        return found;
    }
    bool full() const {
        return true;
    }

    void sort() const {}

    igraph_real_t worstDist() const {
        return radius;
    }

    // Business logic
    // Add the point if it's close enough
    bool addPoint(igraph_real_t dist, igraph_integer_t index) {
        if (dist < radius && index != a && index != b) {
            found += 1;
            if (short_circuit) {
                // Dont continue searching if it only matters if there's one or more.
                return false;
            }
        }
        return true;
    }

};

// Helper result class for listing the points in the intersection of two spheres
// and counting how many are within the lune.
// Used in is_intersection_empty.
// Adapted from code by Szabolcs at https://github.com/szhorvat/IGraphM.
class IntersectionCounts {
    const igraph_real_t radius;                        // Half height of intersection lune
    const igraph_real_t beta_radius;                   // Radius of circles
    const igraph_bool_t short_circuit;                   // Whether to stop after one point has been found
    const igraph_integer_t a, b;                // Edge endpoits; excluded from the count.
    const igraph_vector_t *a_centre, *b_centre; // Circle centres.
    const igraph_matrix_t *points;              // List of points, used to test if points are in range
    size_t count;                               // how many are found.
public:
    // Boilerplate for nanoflann
    using DistanceType = igraph_real_t;
    IntersectionCounts(
        igraph_real_t radius_, igraph_real_t beta_radius_,
        bool short_circuit_,
        igraph_integer_t a, igraph_integer_t b, const igraph_vector_t *a_centre, const igraph_vector_t *b_centre, const igraph_matrix_t *points)
        : radius(radius_), beta_radius(beta_radius_),
          short_circuit(short_circuit_),
          a(a), b(b),
          a_centre(a_centre), b_centre(b_centre), points(points) {
        init();
    }
    void init() {
        clear();
    }

    void clear() {
        count = 0;
    }

    size_t size() const {
        return count;
    }

    bool full() const {
        return true;
    }

    void sort() const {}

    igraph_real_t worstDist() const {
        return radius;
    }
    // Business logic is all contained here.
    // Count the point if it is within the radius from both centres, and if it's not one of the endpoints.
    bool addPoint(igraph_real_t dist, igraph_integer_t index) {
        if (dist < radius && index != a && index != b) {
            igraph_real_t pd1 = vec_ind_sqr_dist(a_centre, index, points);
            igraph_real_t pd2 = vec_ind_sqr_dist(b_centre, index, points);
            if (pd1 < beta_radius && pd2 < beta_radius) {
                count++;
                // Stop searching if it only matters to have at least one point.
                if (short_circuit) {
                    return false;
                }
            }
        }
        return true;
    }
};

// Shrinks type signatures significantly.
template <igraph_integer_t Dimension>
using kdTree = nanoflann::KDTreeSingleIndexAdaptor <
               nanoflann::L2_Adaptor<igraph_real_t, ig_point_adaptor>,
               ig_point_adaptor, Dimension, igraph_integer_t >;

// give a known good superset of edges for a given value of beta.
// In the case of beta < 1, that is a complete graph, for
// beta >= 1 it is the delaunay triangulation of the points.
static igraph_error_t beta_skeleton_edge_superset(igraph_vector_int_t *edges,
        const igraph_matrix_t *points,
        igraph_real_t beta) {

    igraph_integer_t num_points = igraph_matrix_nrow(points);
    igraph_integer_t num_dims   = igraph_matrix_ncol(points);

    if (beta >= 1 && num_points > num_dims) {
        // Large beta and enough points, subset of delaunay.
        IGRAPH_CHECK(igraph_i_delaunay_edges(edges, points));
    } else {
        // Small beta, not subset of delaunay, give complete graph.
        // Or delaunay not calculable due to point count
        // TODO: update when/if delaunay supports small numbers.
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

static inline igraph_real_t calculate_r(igraph_real_t beta) {
    if (beta < 1) {
        return 0.5 / beta;
    }
    return 0.5 * beta;
}


/* -!- centre construction interface -!-
 * igraph_vector_t *a_centre, b_centre : Expects an initialized vector of the
 * correct size, will be written to.
 *
 * igraph_integer_t a, b: the indices of the points,
 *
 * igraph_real_t r: the circles will be constructed with radius r * (distance a-> b)
 *
 * const igraph_matrix_t *points: point set containing the ponits a and b
 */


// construct the centers of the points for lune based beta skeletons with beta
// >= 1. The points lie on the line from a to b, such that the points lie on one
// of the circles.
// formula is a_centre = a + (r-1) * (a - b), similar for b_centre
static igraph_error_t construct_lune_centres(igraph_vector_t *a_centre,
                                      igraph_vector_t *b_centre,
                                      igraph_integer_t a,
                                      igraph_integer_t b,
                                      igraph_real_t r,
                                      const igraph_matrix_t *points) {
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

// construct the centres of the circles for beta < 1, or circle based beta
// skeletons. Since it relies on a 90 degree rotatio, it is only well defined in 2d.
static igraph_error_t construct_perp_centres(igraph_vector_t *a_centre, igraph_vector_t *b_centre, igraph_integer_t a, igraph_integer_t b, igraph_real_t r, const igraph_matrix_t *points) {
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
    // The rotation being x = -y, y = x, 90 degrees counter-clockwise.
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

// Iterate through the edges given, applying the provided filter.
// Currently used with is_intersection_empty and is_union_empty.
// TODO: specialize for multiple dimensions?
template < igraph_error_t filter(igraph_bool_t *result, kdTree < -1 > &tree, igraph_integer_t a, igraph_integer_t b, const igraph_matrix_t *points, igraph_real_t beta) >
static igraph_error_t filter_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points, igraph_real_t beta) {
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

/* -!- Filter interface -!-
 * igraph_bool_t *result        : will be overwritten with the result
 * kdTree <-1> &tree            : used for nanoflann check, should be
 * initialized and filed with points. igraph_integer_t a, b        : the
 * endpoints of the edge being evaluated. igraph_real_t beta           :
 * parameter for beta-skeletons. const igraph_matrix_t *points: matrix
 * containing all the points.
 */

// Test whether the intersection of the two circles as generated
// by centre_positions is empty of points except for a and b.
template <igraph_error_t centre_positions(igraph_vector_t *a_centre,
          igraph_vector_t *b_centre,
          igraph_integer_t a,
          igraph_integer_t b,
          igraph_real_t beta,
          const igraph_matrix_t * points)>
static igraph_error_t is_intersection_empty(
    igraph_bool_t *result,
    kdTree < -1 > &tree,
    igraph_integer_t a,
    igraph_integer_t b,
    const igraph_matrix_t *points,
    igraph_real_t beta) {

    igraph_real_t r = calculate_r(beta);
    igraph_real_t sqr_dist = ind_ind_sqr_distance(a, b, points);

    igraph_vector_t midpoint;
    igraph_integer_t dims = igraph_matrix_ncol(points);


    igraph_vector_t a_centre, b_centre;
    IGRAPH_VECTOR_INIT_FINALLY(&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY(&b_centre, dims);
    IGRAPH_CHECK(centre_positions(&a_centre, &b_centre, a, b, r, points));

    IGRAPH_VECTOR_INIT_FINALLY(&midpoint, dims);


    for (igraph_integer_t i = 0; i < dims; i++) {
        VECTOR(midpoint)[i] =  0.5 * (VECTOR(a_centre)[i] + VECTOR(b_centre)[i]);
    }
    //squared halfheight of lune
    igraph_real_t lune_height = sqr_dist - vec_vec_sqr_dist(&midpoint, &a_centre);

    IntersectionCounts intersections(lune_height, sqr_dist * r * r * (1 + TOLERANCE), true, a, b, &a_centre, &b_centre, points);

    tree.findNeighbors(intersections, VECTOR(midpoint));

    *result = intersections.size() == 0;

    igraph_vector_destroy(&midpoint);
    igraph_vector_destroy(&a_centre);
    igraph_vector_destroy(&b_centre);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

// Test whether the union of the circles given by centre_positions
//  is empty of other points except for a and b.
template <igraph_error_t centre_positions(igraph_vector_t *a_centre,
          igraph_vector_t *b_centre,
          igraph_integer_t a,
          igraph_integer_t b,
          igraph_real_t beta,
          const igraph_matrix_t * points)>
static igraph_error_t is_union_empty(
    igraph_bool_t *result,
    kdTree < -1 > &tree,
    igraph_integer_t a,
    igraph_integer_t b,
    const igraph_matrix_t *points,
    igraph_real_t beta) {

    igraph_integer_t dims = igraph_matrix_ncol(points);

    igraph_real_t sqr_dist = ind_ind_sqr_distance(a, b, points);
    igraph_real_t r = calculate_r(beta);
    igraph_vector_t a_centre, b_centre;



    IGRAPH_VECTOR_INIT_FINALLY(&a_centre, dims);
    IGRAPH_VECTOR_INIT_FINALLY(&b_centre, dims);
    IGRAPH_CHECK(centre_positions(&a_centre, &b_centre, a, b, r, points));

    NeighborCounts neighbor_search(sqr_dist * r * r * (1 + TOLERANCE), a, b, true);

    tree.findNeighbors(neighbor_search, VECTOR(a_centre));


    if (neighbor_search.size() > 0) {
        *result = false;
    } else {
        neighbor_search.clear();
        tree.findNeighbors(neighbor_search, VECTOR(b_centre));
        *result = neighbor_search.size() == 0;
    }
    igraph_vector_destroy(&a_centre);
    igraph_vector_destroy(&b_centre);
    IGRAPH_FINALLY_CLEAN(2);
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
    // determine filter required based on beta.
    if (beta >= 1) {
        IGRAPH_CHECK(filter_edges<is_intersection_empty<construct_lune_centres>>(&potential_edges, points, beta));
    } else {
        if (igraph_matrix_ncol(points) != 2) {
            IGRAPH_ERROR("Beta skeletons with beta < 1 are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
        }

        IGRAPH_CHECK(filter_edges<is_intersection_empty<construct_perp_centres>>(&potential_edges, points, beta));
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
        IGRAPH_CHECK(filter_edges<is_union_empty<construct_perp_centres>>(&potential_edges, points, beta));
    } else {
        IGRAPH_CHECK(filter_edges<is_intersection_empty<construct_perp_centres>>(&potential_edges, points, beta));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

// Derived from code by Szabolcs at github.com/SzHorvat/IGraphM

class BetaFinder {
    const igraph_real_t max_beta;
    const igraph_real_t tol;
    const igraph_integer_t ai, bi;
    const igraph_matrix_t *ps;
    const igraph_real_t ab2;

    igraph_real_t smallest_beta;
    igraph_real_t max_radius;

public:
    using DistanceType = igraph_real_t;
    using IndexType = igraph_integer_t;
    BetaFinder(igraph_real_t max_beta, igraph_real_t tol, igraph_integer_t v1, igraph_integer_t v2, const igraph_matrix_t *ps) :
        max_beta(max_beta), tol(tol), ai(v1), bi(v2), ps(ps),
        ab2(ind_ind_sqr_distance(ai, bi, ps)) {
        init();
    }

    void init() {
        clear();
    }

    void clear() {
        smallest_beta = IGRAPH_INFINITY;
        max_radius = luneHalfHeight2(max_beta);
    }

    bool full() const {
        return true;
    }

    size_t size() const {
        return 1;
    }

    igraph_real_t luneHalfHeight2(igraph_real_t beta) const {
        if (beta == 0) {
            return 0;
        }
        return (ab2 / 4) * (2 * beta - 1);
    }
    // Calculate the beta at which point index would make the edge a-b diappear.
    // If calculated beta is under 1, it would not appear in the gabriel graph, so a 0 is returned
    // which is to be interpreted as the edge being missing.
    igraph_real_t pointBeta(igraph_integer_t index) const {
        igraph_real_t ap2 = ind_ind_sqr_distance(ai, index, ps);
        igraph_real_t bp2 = ind_ind_sqr_distance(bi, index, ps);

        if (ap2 > bp2) {
            std::swap(ap2, bp2);
        }

        igraph_real_t denom = ab2 + ap2 - bp2;

        if (denom <= 0) {
            return IGRAPH_INFINITY;
        }

        igraph_real_t beta = 2 * ap2 / denom;

        if (beta < 1 + tol) {
            printf("XXX Dropping edge (%d, %d) with index=%d, beta=%g, denom=%g, ab2=%g, ap2=%g, bp2=%g\n", (int) ai, (int) bi, (int) index, beta, denom, ab2, ap2, bp2);
        }

        return beta < 1 + tol ? 0 : beta;
    }

    bool addPoint(igraph_real_t dist, igraph_integer_t index) {

        if (dist < max_radius) {
            igraph_real_t beta = pointBeta(index);
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
    igraph_real_t thresholdBeta() const {
        return smallest_beta;
    }
};

/**
 * \function igraph_beta_weighted_gabriel_graph
 *
 *  \brief Computes a weighted graph where edge weight represent the beta value at which the edge would disappear.
 *
 * \experimental
 *
 * This function computes a gabriel graph, where the weight of each edge represents
 *   the value of beta at which the edge would disappear in a lune based beta skeleton.
 * Larger values of max_beta are slower to calculate.
 *
 * \param graph A pointer to the graph that will be created.
 * \param edge_weights Will contain the edge weights corresponding to the edge indices from the graph.
 * \param points A matrix containing the points that will be used.
 *     Each row is a point, dimensionality is inferred from the column count.
 *    There must be no duplicate points.
 * \param max_beta Maximum value of beta to search to, higher values will be represented as infinity.
 * \return Error code.
 *
 * \sa \ref igraph_lune_beta_skeleton() or \ref igraph_circle_beta_skeleton() to generate a graph with a given value of beta.
 *
 * Time complexity: TODO
 */
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
        BetaFinder finder(max_beta, TOLERANCE, VECTOR(edges)[2 * i], VECTOR(edges)[2 * i + 1], points);

        for (igraph_integer_t axis = 0; axis < dim; axis++) {
            VECTOR(midpoint)[axis] = 0.5 * (
                                         MATRIX(*points, VECTOR(edges)[2 * i],     axis)
                                         + MATRIX(*points, VECTOR(edges)[2 * i + 1], axis)
                                     );
        }

        tree.findNeighbors(finder, VECTOR(midpoint));
        VECTOR(*edge_weights)[i] = finder.thresholdBeta();

    }

    // Filter out edges with beta == 0.
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
