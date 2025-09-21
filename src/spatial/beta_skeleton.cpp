/*
   igraph library.
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

#include "igraph_spatial.h"

#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "core/exceptions.h"
#include "spatial/nanoflann_internal.hpp"
#include "spatial/spatial_internal.h"
#include <cfloat>

#define TOLERANCE (128 * DBL_EPSILON)

// Some methods to get distances between two vectors, including when one or both are embedded in a matrix.
static inline igraph_real_t ind_ind_sqr_distance(igraph_int_t a, igraph_int_t b,
        const igraph_matrix_t *points) {
    igraph_real_t distance = 0;
    igraph_int_t dims = igraph_matrix_ncol(points);
    for (igraph_int_t i = 0; i < dims; i++) {
        igraph_real_t temp = MATRIX(*points, a, i) - MATRIX(*points, b, i);
        distance += temp * temp;
    }
    return distance;
}

static inline igraph_real_t vec_vec_sqr_dist(const igraph_vector_t *a, const igraph_vector_t *b) {
    igraph_real_t distance = 0;
    igraph_int_t dims = igraph_vector_size(a);
    for (igraph_int_t i = 0; i < dims; i++) {
        igraph_real_t temp = VECTOR(*a)[i] - VECTOR(*b)[i];
        distance += temp * temp;
    }
    return distance;
}

static inline igraph_real_t vec_vec_sqr_dist(const std::vector<igraph_real_t> &a, const std::vector<igraph_real_t> &b) {
    igraph_real_t distance = 0;
    igraph_int_t dims = a.size();
    for (igraph_int_t i = 0; i < dims; i++) {
        igraph_real_t temp = a[i] - b[i];
        distance += temp * temp;
    }
    return distance;
}

static inline igraph_real_t vec_ind_sqr_dist(const igraph_vector_t *a, const igraph_int_t b, const igraph_matrix_t *points) {
    igraph_real_t distance = 0;
    igraph_int_t dims = igraph_matrix_ncol(points);
    for (igraph_int_t i = 0; i < dims; i++) {
        igraph_real_t temp = VECTOR(*a)[i] - MATRIX(*points, b, i);
        distance += temp * temp;
    }
    return distance;
}
static inline igraph_real_t vec_ind_sqr_dist(const std::vector<igraph_real_t> &a, const igraph_int_t b, const igraph_matrix_t *points) {
    igraph_real_t distance = 0;
    igraph_int_t dims = igraph_matrix_ncol(points);
    for (igraph_int_t i = 0; i < dims; i++) {
        igraph_real_t temp = a[i] - MATRIX(*points, b, i);
        distance += temp * temp;
    }
    return distance;
}

// Adapted from code by Szabolcs at https://github.com/szhorvat/IGraphM
// Used in is_union_empty
class NeighborCounts {
    const igraph_real_t radius; // L2 search radius
    const igraph_int_t a, b; // edge endpoints; excluded from the count
    igraph_bool_t short_circuit;
    igraph_int_t found;

public:
    // Boilerplate for nanoflann
    using DistanceType = igraph_real_t;
    NeighborCounts(igraph_real_t radius, igraph_int_t a, igraph_int_t b, igraph_bool_t short_circuit)
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
    bool addPoint(igraph_real_t dist, igraph_int_t index) {
        if (dist < radius && index != a && index != b) {
            found += 1;
            if (short_circuit) {
                // Don't continue searching if it only matters if there's one or more.
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
    const igraph_real_t radius;                            // Half height of intersection lune
    const igraph_real_t beta_radius;                       // Radius of circles
    const igraph_bool_t short_circuit;                     // Whether to stop after one point has been found
    const igraph_int_t a, b;                               // Edge endpoints; excluded from the count.
    const std::vector<igraph_real_t> &a_center, &b_center; // Circle centers.
    const igraph_matrix_t *points;                         // List of points, used to test if points are in range
    size_t count;                                          // how many are found.
public:
    // Boilerplate for nanoflann
    using DistanceType = igraph_real_t;
    IntersectionCounts(
        igraph_real_t radius_, igraph_real_t beta_radius_,
        bool short_circuit_,
        igraph_int_t a, igraph_int_t b, const std::vector<igraph_real_t> &a_center, const std::vector<igraph_real_t> &b_center, const igraph_matrix_t *points)
        : radius(radius_), beta_radius(beta_radius_),
          short_circuit(short_circuit_),
          a(a), b(b),
          a_center(a_center), b_center(b_center), points(points) {
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
    // Count the point if it is within the radius from both centers, and if it's not one of the endpoints.
    bool addPoint(igraph_real_t dist, igraph_int_t index) {
        if (dist < radius && index != a && index != b) {
            igraph_real_t pd1 = vec_ind_sqr_dist(a_center, index, points);
            igraph_real_t pd2 = vec_ind_sqr_dist(b_center, index, points);
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
template<igraph_int_t Dimension>
using KDTree = nanoflann::KDTreeSingleIndexAdaptor <
               nanoflann::L2_Adaptor<igraph_real_t, ig_point_adaptor>,
               ig_point_adaptor, Dimension, igraph_int_t >;

// give a known good superset of edges for a given value of beta.
// In the case of beta < 1, that is a complete graph, for
// beta >= 1 it is the delaunay triangulation of the points.
static igraph_error_t beta_skeleton_edge_superset(igraph_vector_int_t *edges,
        const igraph_matrix_t *points,
        igraph_real_t beta) {

    igraph_int_t num_points = igraph_matrix_nrow(points);
    igraph_int_t num_dims   = igraph_matrix_ncol(points);

    if (beta >= 1 && num_points > num_dims) {
        // Large beta and enough points, subset of delaunay.
        IGRAPH_CHECK(igraph_i_delaunay_edges(edges, points));
    } else {
        // Small beta, not subset of Delaunay, give complete graph.
        // Or Delaunay not calculable due to point count
        // TODO: update when/if Delaunay supports small numbers.
        igraph_int_t numpoints = igraph_matrix_nrow(points);
        for (igraph_int_t a = 0; a < numpoints - 1; a++) {
            for (igraph_int_t b = a + 1; b < numpoints; b++) {
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


/* -!- center construction interface -!-
 * igraph_vector_t *a_center, b_center : Expects an initialized vector of the
 * correct size, will be written to.
 *
 * igraph_int_t a, b: the indices of the points,
 *
 * igraph_real_t r: the circles will be constructed with radius r * (distance a-> b)
 *
 * const igraph_matrix_t *points: point set containing the points a and b
 */

typedef void CenterConstructor(
    std::vector<igraph_real_t> &a_centre,
    std::vector<igraph_real_t> &b_centre,
    igraph_int_t a,
    igraph_int_t b,
    igraph_real_t beta,
    const igraph_matrix_t *points);

// construct the centers of the points for lune based beta skeletons with beta
// >= 1. The points lie on the line from a to b, such that the points lie on one
// of the circles.
// formula is a_center = a + (r-1) * (a - b), similar for b_center
static void construct_lune_centers(std::vector<igraph_real_t> &a_centre,
                                   std::vector<igraph_real_t> &b_centre,
                                   igraph_int_t a,
                                   igraph_int_t b,
                                   igraph_real_t r,
                                   const igraph_matrix_t *points) {
    igraph_int_t dims = igraph_matrix_ncol(points);
    a_centre.resize(dims);
    b_centre.resize(dims);
    for (igraph_int_t i = 0; i < dims; i++) {
        a_centre[i] = MATRIX(*points, a, i) + (r - 1) * (MATRIX(*points, a, i) - MATRIX(*points, b, i));
        b_centre[i] = MATRIX(*points, b, i) + (r - 1) * (MATRIX(*points, b, i) - MATRIX(*points, a, i));
    }
}

// construct the centers of the circles for beta < 1, or circle based beta
// skeletons.
// Since it relies on a 90-degree rotation around the axis perpendicular
// to the line AB, it is only well-defined in 2d.
static void construct_perp_centers(std::vector<igraph_real_t> &a_centre,
                                   std::vector<igraph_real_t> &b_centre,
                                   igraph_int_t a,
                                   igraph_int_t b,
                                   igraph_real_t r,
                                   const igraph_matrix_t *points) {
    igraph_real_t mid[2], perp[2];

    for (igraph_int_t i = 0; i < 2; i++) {
        mid[i]  = (MATRIX(*points, a, i) + MATRIX(*points, b, i)) * 0.5;
        perp[i] = (MATRIX(*points, a, i) - MATRIX(*points, b, i)) * sqrt(r * r - 0.25);
    }

    // Since this is only well-defined for 2d, a manual 90-degree rotation works and is simpler.
    // The rotation being x = -y, y = x, 90 degrees counter-clockwise.
    igraph_real_t temp = perp[0];
    perp[0] = - perp[1];
    perp[1] = temp;

    a_centre.resize(2);
    b_centre.resize(2);
    for (igraph_int_t i = 0; i < 2; i++) {
        a_centre[i] = mid[i] + perp[i];
        b_centre[i] = mid[i] - perp[i];
    }
}

/* -!- Filter interface -!-
 * igraph_bool_t *result        : will be overwritten with the result
 * kdTree <-1> &tree            : used for nanoflann check, should be
 * initialized and filed with points. igraph_int_t a, b        : the
 * endpoints of the edge being evaluated. igraph_real_t beta           :
 * parameter for beta-skeletons. const igraph_matrix_t *points: matrix
 * containing all the points.
 */
typedef bool FilterFunc(
    const KDTree < -1 > &tree,
    igraph_int_t a,
    igraph_int_t b,
    const igraph_matrix_t *points,
    igraph_real_t beta
);


// Test whether the intersection of the two circles as generated
// by center_positions is empty of points except for a and b.
template<CenterConstructor center_positions, igraph_bool_t is_closed>
static bool is_intersection_empty(
    const KDTree < -1 > &tree,
    igraph_int_t a,
    igraph_int_t b,
    const igraph_matrix_t *points,
    igraph_real_t beta) {

    igraph_real_t r = calculate_r(beta);
    igraph_real_t sqr_dist = ind_ind_sqr_distance(a, b, points);

    std::vector<igraph_real_t> midpoint;
    igraph_int_t dims = igraph_matrix_ncol(points);

    std::vector<igraph_real_t> a_centre, b_centre;
    center_positions(a_centre, b_centre, a, b, r, points);

    midpoint.resize(dims);
    for (igraph_int_t i = 0; i < dims; i++) {
        midpoint[i] =  0.5 * (a_centre[i] + b_centre[i]);
    }
    // squared half-height of lune
    igraph_real_t lune_height = sqr_dist - vec_vec_sqr_dist(midpoint, a_centre);
    igraph_real_t tol = is_closed ? 1 + TOLERANCE : 1 - TOLERANCE;
    IntersectionCounts intersections(lune_height * tol * tol, sqr_dist * r * r * tol * tol, true, a, b, a_centre, b_centre, points);

    tree.findNeighbors(intersections, midpoint.data());
    return intersections.size() == 0;
}

// Test whether the union of the circles given by center_positions
//  is empty of other points except for a and b.
template<CenterConstructor center_positions>
static bool is_union_empty(
    const KDTree < -1 > &tree,
    igraph_int_t a,
    igraph_int_t b,
    const igraph_matrix_t *points,
    igraph_real_t beta) {

    igraph_real_t sqr_dist = ind_ind_sqr_distance(a, b, points);
    igraph_real_t r = calculate_r(beta);
    std::vector<igraph_real_t> a_centre, b_centre;

    center_positions(a_centre, b_centre, a, b, r, points);

    NeighborCounts neighbor_search(sqr_dist * r * r * (1 + TOLERANCE), a, b, true);

    tree.findNeighbors(neighbor_search, a_centre.data());

    if (neighbor_search.size() > 0) {
        return false;
    } else {
        neighbor_search.clear();
        tree.findNeighbors(neighbor_search, b_centre.data());
        return neighbor_search.size() == 0;
    }
}

// Iterate through the edges given, applying the provided filter.
// Currently used with is_intersection_empty and is_union_empty.
// TODO: specialize for multiple dimensions?
template<FilterFunc filter>
static igraph_error_t filter_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_int_t available_edges = igraph_vector_int_size(edges);
    igraph_int_t added_edges = 0;
    ig_point_adaptor adaptor(points);
    igraph_int_t dim = igraph_matrix_ncol(points);
    KDTree < -1 > tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();
    for (igraph_int_t i = 0; i * 2  < available_edges; i++) {
        if (filter(tree, VECTOR(*edges)[2 * i], VECTOR(*edges)[2 * i + 1], points, beta)) {
            VECTOR(*edges)[added_edges * 2]     = VECTOR(*edges)[i * 2];
            VECTOR(*edges)[added_edges * 2 + 1] = VECTOR(*edges)[i * 2 + 1];
            added_edges += 1;
        }
    }
    IGRAPH_CHECK(igraph_vector_int_resize(edges, added_edges * 2));
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_lune_beta_skeleton
 * \brief The lune based β-skeleton of a spatial point set.
 *
 * \experimental
 *
 * This function constructs the lune-based β-skeleton of an n-dimensional
 * spatial point set.
 *
 * </para><para>
 * A larger β results in a larger region, and a sparser graph.
 * Values of β &lt; 1 are only supported in 2D, and are considerably slower.
 *
 * </para><para>
 * The Gabriel graph is a special case of beta skeleton where <code>β = 1</code>.
 *
 * </para><para>
 * The Relative Neighborhood graph is a special case of beta skeleton where
 * β approaches
 *
 * \param graph A pointer to the graph that will be created.
 * \param points A matrix containing the points that will be used to create the
 *     graph. Each row is a point, dimensionality is inferred from the column count.
 *
 * \return Error code.
 *
 * Time Complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_lune_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    // determine filter required based on beta.
    if (beta >= 1) {
        IGRAPH_CHECK((filter_edges<is_intersection_empty<construct_lune_centers, true>>(&potential_edges, points, beta)));
    } else {
        if (igraph_matrix_ncol(points) != 2) {
            IGRAPH_ERROR("Beta skeletons with beta < 1 are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
        }

        IGRAPH_CHECK((filter_edges<is_intersection_empty<construct_perp_centers, true>>(&potential_edges, points, beta)));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_circle_beta_skeleton
 * \brief The circle based β-skeleton of a 2D spatial point set.
 *
 * \experimental
 *
 * This function constructs the circle based β-skeleton of a 2D spatial point set.
 *
 * </para><para>
 * A larger \p beta value results in a larger region, and a sparser graph.
 * Values of beta &lt; 1 are considerably slower
 *
 * \param graph A pointer to the graph that will be created.
 * \param points An n-by-2 matrix containing the points that will be used to
 *    create the graph. Each row is a point.
 * \param beta A positive real value used to parameterize the graph.
 * \return Error code.
 *
 * Time Complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_circle_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    if (igraph_matrix_ncol(points) != 2) {
        IGRAPH_ERROR("Circle based beta skeletons are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
    }

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    if (beta >= 1) {
        IGRAPH_CHECK(filter_edges<is_union_empty<construct_perp_centers>>(&potential_edges, points, beta));
    } else {
        IGRAPH_CHECK((filter_edges<is_intersection_empty<construct_perp_centers, true>>(&potential_edges, points, beta)));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

// Derived from code by Szabolcs at github.com/szhorvat/IGraphM

class BetaFinder {
    const igraph_real_t max_beta;
    const igraph_real_t tol;
    const igraph_int_t ai, bi;
    const igraph_matrix_t *ps;
    const igraph_real_t ab2;

    igraph_real_t smallest_beta;
    igraph_real_t max_radius;

public:
    using DistanceType = igraph_real_t;
    using IndexType = igraph_int_t;
    BetaFinder(igraph_real_t max_beta, igraph_real_t tol, igraph_int_t v1, igraph_int_t v2, const igraph_matrix_t *ps) :
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

    // Calculate the beta at which point index would make the edge a-b disappear.
    // If calculated beta is under 1, it would not appear in the gabriel graph, so a 0 is returned
    // which is to be interpreted as the edge being missing.
    igraph_real_t pointBeta(igraph_int_t index) const {
        if (index == ai || index == bi) {
            return IGRAPH_INFINITY;
        }

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
        return beta < 1 + tol ? 0 : beta;
    }

    bool addPoint(igraph_real_t dist, igraph_int_t index) {
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
 * \brief A Gabriel graph, with edges weighted by the β value at which it disappears.
 *
 * \experimental
 *
 * This function generates a Gabriel graph, and for each edge of this graph it
 * computes the threshold β value at which the edge ceases to be part of the
 * lune-based β-skeleton. For edges that continue to be part of β-skeletons
 * for arbitrarily large β, \c IGRAPH_INFINITΥ is returned.
 *
 * </para><para>
 * The \p max_beta cutoff parameter controls the largest β value to consider
 * For edges that persist above this β value, \c IGRAPH_INFINITΥ is returned.
 * This parameter serves to improve performance: the smaller this cutoff,
 * the faster the computation. Pass \c IGRAPH_INFINITY to use no cutoff.
 *
 * \param graph A pointer to the graph that will be created.
 * \param weights Will contain the edge weights corresponding to the edge
 *    indices from the graph.
 * \param points A matrix containing the points that will be used.
 *    Each row is a point, dimensionality is inferred from the column count.
 *    There must be no duplicate points.
 * \param max_beta Maximum value of beta to search to, higher values will be
 *    represented as \c IGRAPH_INFINITY.
 * \return Error code.
 *
 * \sa \ref igraph_lune_beta_skeleton() or \ref igraph_circle_beta_skeleton()
 * to generate a graph with a given value of beta; \ref igraph_gabriel_graph()
 * to only generate a Gabriel graph, without edge weights.
 *
 *
 * Time Complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set. Though large values of max_beta
 * can cause long run times if there are edges that disappear only at large betas.
 *
 */
igraph_error_t igraph_beta_weighted_gabriel_graph(
    igraph_t *graph,
    igraph_vector_t *weights,
    const igraph_matrix_t *points,
    igraph_real_t max_beta) {

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    igraph_int_t dim = igraph_matrix_ncol(points);
    igraph_int_t point_count = igraph_matrix_nrow(points);
    ig_point_adaptor adaptor(points);

    IGRAPH_CHECK(igraph_i_delaunay_edges(&edges, points));
    igraph_int_t edge_count = igraph_vector_int_size(&edges) / 2;

    IGRAPH_CHECK(igraph_vector_resize(weights, edge_count));

    KDTree < -1 > tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

    igraph_vector_t midpoint;
    IGRAPH_VECTOR_INIT_FINALLY(&midpoint, dim);

    for (igraph_int_t i = 0; i < edge_count; i++) {
        BetaFinder finder(max_beta, TOLERANCE, VECTOR(edges)[2 * i], VECTOR(edges)[2 * i + 1], points);

        for (igraph_int_t axis = 0; axis < dim; axis++) {
            VECTOR(midpoint)[axis] = 0.5 * (
                                         MATRIX(*points, VECTOR(edges)[2 * i],     axis)
                                         + MATRIX(*points, VECTOR(edges)[2 * i + 1], axis)
                                     );
        }

        tree.findNeighbors(finder, VECTOR(midpoint));
        VECTOR(*weights)[i] = finder.thresholdBeta();
    }

    // Filter out edges with beta == 0.
    igraph_int_t added_edges = 0;
    for (igraph_int_t i = 0; i < edge_count; i++) {
        if (VECTOR(*weights)[i] != 0) {
            VECTOR(*weights)[added_edges] = VECTOR(*weights)[i];
            VECTOR(edges)[added_edges * 2] = VECTOR(edges)[i * 2];
            VECTOR(edges)[added_edges * 2 + 1] = VECTOR(edges)[i * 2 + 1];
            added_edges += 1;
        }
    }

    IGRAPH_CHECK(igraph_vector_int_resize(&edges, added_edges * 2));
    IGRAPH_CHECK(igraph_vector_resize(weights, added_edges));

    IGRAPH_CHECK(igraph_create(graph, &edges, point_count, false));

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&midpoint);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_gabriel_graph
 * \brief The Gabriel graph of a point set.
 *
 * \experimental
 *
 * In the Gabriel graph of a point set, two points A and B are connected if
 * there is no other point C within the closed ball of which AB is a diameter.
 * The Gabriel graph is connected, and in 2D it is planar. igraph supports
 * computing the Gabriel graph of arbitrary dimensional point sets.
 *
 * </para><para>
 * The Gabriel graph is a special case of lune-based and circle-based β-skeletons
 * with <code>β=1</code>.
 *
 * \param graph A pointer to the graph to be created.
 * \param points The point set that will be used. Each row is a point,
 *    dimensionality is inferred from column count.
 *
 * \return Error Code.
 * \sa The Gabriel graph is a special case of
 *     \ref igraph_lune_beta_skeleton() and \ref igraph_circle_beta_skeleton()
 *     where <code>β = 1</code>.
 *
 * Time Complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_gabriel_graph(igraph_t *graph, const igraph_matrix_t *points) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, 1));

    IGRAPH_CHECK((filter_edges<is_intersection_empty<construct_lune_centers, true>>(&potential_edges, points, 1)));

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_relative_neighborhood_graph
 * \brief The relative neighborhood graph of a point set.
 *
 * \experimental
 *
 * The relative neighborhood graph is constructed from a set of points in space.
 * Two points A and B are connected if and only if there is no other point C so
 * that AC &lt; AB and BC &lt; AB, with the inequalities being strict.
 *
 * </para><para>
 * Most authors define the relative neighborhood graph to coincide with a
 * lune-based β-skeleton for <code>β = 2</code>. In igraph, there is a subtle
 * difference: the <code>β = 2</code> skeleton connects points A and B when there
 * is no point C so that AC &lt;= AB and BC &lt;= AB. Therefore, three points
 * forming an equilateral triangle are connected in the relative neighborhood graph,
 * but disconnected in the <code>β = 2</code> skeleton.
 *
 * </para><para>
 * With these definitions, the relative neighborhood graph is always connected,
 * while the <code>β = 2</code> skeleton is always triangle-free.
 *
 * \param graph A pointer to the graph that will be created.
 * \param points The point set that will be used, each row is a point.
 *     Dimensionality is inferred from the column number.
 * \return Error code.
 *
 * \sa \ref igraph_lune_beta_skeleton() to compute the lune based β-skeleton
 * for <code>β = 2</code> or other β values.
 *
 * Time Complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_relative_neighborhood_graph(igraph_t *graph, const igraph_matrix_t *points) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, 1));

    IGRAPH_CHECK((filter_edges<is_intersection_empty<construct_lune_centers, false>>(&potential_edges, points, 2)));

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}
