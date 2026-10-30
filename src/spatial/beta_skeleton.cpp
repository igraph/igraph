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
#include "core/interruption.h"
#include "spatial/nanoflann_internal.hpp"
#include "spatial/spatial_internal.h"
#include <algorithm>
#include <cfloat>
#include <cmath>

#define TOLERANCE (128 * DBL_EPSILON)

// The exclusion regions of beta skeletons are bounded by spheres through the
// endpoints a and b of the tested edge. Membership is tested using the
// differences p - a, p - b and b - a, instead of constructing the sphere centres
// in absolute coordinates. Differences of nearby points have small relative
// errors even when the points are far from the origin, so points on the boundary
// of a closed region are detected reliably. The k-d tree is only used to find
// candidate points, so its search radii include a margin for rounding errors.

// Squared radius of a k-d tree search ball around 'centre' that contains all
// points within sqrt(sqr_radius) of the exact centre. The margin accounts for
// rounding errors in computing the centre from the input coordinates, which are
// proportional to the magnitude of the coordinates, and in nanoflann's distances.
static igraph_real_t search_sqr_radius(const std::vector<igraph_real_t> &centre, igraph_real_t sqr_radius) {
    igraph_real_t magnitude = 0;
    for (const igraph_real_t c : centre) {
        magnitude += fabs(c);
    }
    const igraph_real_t radius = sqrt(sqr_radius) * (1 + TOLERANCE) + 16 * DBL_EPSILON * magnitude;
    return radius * radius;
}

// The allowed value of |p - c|^2 - R^2 for points p in a ball of squared radius
// sqr_radius = R^2. It is positive for closed regions and negative for open ones,
// and corresponds to scaling the radius by 1 + TOLERANCE or 1 - TOLERANCE.
static inline igraph_real_t boundary_band(igraph_real_t sqr_radius, bool is_closed) {
    const igraph_real_t tol = is_closed ? 1 + TOLERANCE : 1 - TOLERANCE;
    return (tol * tol - 1) * sqr_radius;
}

static inline igraph_real_t calculate_r(igraph_real_t beta) {
    if (beta < 1) {
        return 0.5 / beta;
    }
    return 0.5 * beta;
}

// Data common to the exclusion regions of an edge a-b. Objects are reused
// for all edges to avoid memory allocations; call set_edge() for each edge.
class EdgeRegion {
protected:
    const igraph_matrix_t *points;
    const igraph_int_t dim;
    igraph_int_t a = -1, b = -1;
    std::vector<igraph_real_t> u; // b - a
    igraph_real_t uu = 0;         // |b - a|^2

    explicit EdgeRegion(const igraph_matrix_t *points) :
        points(points), dim(igraph_matrix_ncol(points)), u(dim) { }

    void set_edge(igraph_int_t a_, igraph_int_t b_) {
        a = a_;
        b = b_;
        uu = 0;
        for (igraph_int_t i = 0; i < dim; i++) {
            u[i] = MATRIX(*points, b, i) - MATRIX(*points, a, i);
            uu += u[i] * u[i];
        }
    }

    void midpoint(std::vector<igraph_real_t> &centre) const {
        centre.resize(dim);
        for (igraph_int_t i = 0; i < dim; i++) {
            centre[i] = 0.5 * (MATRIX(*points, a, i) + MATRIX(*points, b, i));
        }
    }

    // Computes |p-a|^2, |p-b|^2 and (p-a).(p-b).
    void local_products(igraph_int_t p, igraph_real_t &vv, igraph_real_t &ww, igraph_real_t &vw) const {
        vv = ww = vw = 0;
        for (igraph_int_t i = 0; i < dim; i++) {
            const igraph_real_t v = MATRIX(*points, p, i) - MATRIX(*points, a, i);
            const igraph_real_t w = MATRIX(*points, p, i) - MATRIX(*points, b, i);
            vv += v * v;
            ww += w * w;
            vw += v * w;
        }
    }
};

// The lune of a-b for beta >= 1, in any dimension: the intersection of two
// balls of radius beta/2 |b - a| centred on the line through a and b, so that
// b is on the boundary of the first ball and a on the boundary of the second.
// With v = p - a and w = p - b, p is in the first ball when
// beta v.w - (beta - 1) |w|^2 <= 0, and in the second when
// beta v.w - (beta - 1) |v|^2 <= 0. These expressions equal |p - c|^2 - R^2
// for the ball centre c and radius R. For beta = 1, both reduce to v.w <= 0,
// the test for the Gabriel ball.
template<bool is_closed>
class LuneRegion : public EdgeRegion {
    const igraph_real_t beta;
    igraph_real_t band = 0;

public:
    static const int query_count = 1;

    LuneRegion(const igraph_matrix_t *points, igraph_real_t beta) :
        EdgeRegion(points), beta(beta) { }

    void set_edge(igraph_int_t a, igraph_int_t b) {
        EdgeRegion::set_edge(a, b);
        band = boundary_band(0.25 * beta * beta * uu, is_closed);
    }

    bool contains(igraph_int_t p) const {
        if (p == a || p == b) {
            return false;
        }
        igraph_real_t vv, ww, vw;
        local_products(p, vv, ww, vw);
        const igraph_real_t f1 = beta * vw - (beta - 1) * ww, f2 = beta * vw - (beta - 1) * vv;
        return is_closed ? (f1 <= band && f2 <= band) : (f1 < band && f2 < band);
    }

    // The lune lies within sqrt(2 beta - 1) |b - a| / 2 of the midpoint.
    igraph_real_t query(int, std::vector<igraph_real_t> &centre) const {
        midpoint(centre);
        return 0.25 * (2 * beta - 1) * uu + std::max(band, 0.0);
    }
};

// Regions bounded by the two circles through a and b with radius R = r |b - a|.
// Their centres are m + s J(b-a) and m - s J(b-a), where m is the midpoint of a
// and b, J is a rotation by 90 degrees and s = sqrt(r^2 - 1/4). For these centres,
// |p - c|^2 - R^2 equals (p-a).(p-b) - 2 s cross(b-a, p-a) and
// (p-a).(p-b) + 2 s cross(b-a, p-a). Only defined in 2D.
class CircleRegion : public EdgeRegion {
public:
    void set_edge(igraph_int_t a, igraph_int_t b) {
        EdgeRegion::set_edge(a, b);
        band = boundary_band(r * r * uu, true);
    }

protected:
    const igraph_real_t r, s;
    igraph_real_t band = 0;

    CircleRegion(const igraph_matrix_t *points, igraph_real_t beta) :
        EdgeRegion(points), r(calculate_r(beta)), s(sqrt(r * r - 0.25)) { }

    // Computes (p-a).(p-b) and cross(b-a, p-a).
    void circle_products(igraph_int_t p, igraph_real_t &dot, igraph_real_t &cross) const {
        const igraph_real_t v0 = MATRIX(*points, p, 0) - MATRIX(*points, a, 0);
        const igraph_real_t v1 = MATRIX(*points, p, 1) - MATRIX(*points, a, 1);
        const igraph_real_t w0 = MATRIX(*points, p, 0) - MATRIX(*points, b, 0);
        const igraph_real_t w1 = MATRIX(*points, p, 1) - MATRIX(*points, b, 1);
        dot = v0 * w0;
        dot += v1 * w1;
        cross = u[0] * v1 - u[1] * v0;
    }
};

// The closed intersection of the two discs, used for both lune-based and
// circle-based skeletons with beta < 1.
class CircleIntersectionRegion : public CircleRegion {
public:
    static const int query_count = 1;

    CircleIntersectionRegion(const igraph_matrix_t *points, igraph_real_t beta) :
        CircleRegion(points, beta) { }

    bool contains(igraph_int_t p) const {
        if (p == a || p == b) {
            return false;
        }
        igraph_real_t dot, cross;
        circle_products(p, dot, cross);
        return dot + 2 * s * fabs(cross) <= band;
    }

    // The intersection lies within |b - a| / 2 of the midpoint.
    igraph_real_t query(int, std::vector<igraph_real_t> &centre) const {
        midpoint(centre);
        return 0.25 * uu + band;
    }
};

// The closed union of the two discs, used for circle-based skeletons with beta >= 1.
class CircleUnionRegion : public CircleRegion {
public:
    static const int query_count = 2;

    CircleUnionRegion(const igraph_matrix_t *points, igraph_real_t beta) :
        CircleRegion(points, beta) { }

    bool contains(igraph_int_t p) const {
        if (p == a || p == b) {
            return false;
        }
        igraph_real_t dot, cross;
        circle_products(p, dot, cross);
        return dot - 2 * s * fabs(cross) <= band;
    }

    // Query k searches the disc centred at m + s J(b-a) or m - s J(b-a).
    igraph_real_t query(int k, std::vector<igraph_real_t> &centre) const {
        midpoint(centre);
        const igraph_real_t sign = k == 0 ? 1 : -1;
        centre[0] -= sign * s * u[1];
        centre[1] += sign * s * u[0];
        return r * r * uu + band;
    }
};

// nanoflann result set that stops at the first point in the region.
template<typename Region>
class RegionSearch {
    const Region &region;
    const igraph_real_t sqr_radius;
    bool found;

public:
    // Boilerplate for nanoflann
    using DistanceType = igraph_real_t;
    using IndexType = igraph_int_t;

    RegionSearch(const Region &region, igraph_real_t sqr_radius) :
        region(region), sqr_radius(sqr_radius) {
        init();
    }

    void init() {
        clear();
    }

    void clear() {
        found = false;
    }

    size_t size() const {
        return found;
    }

    bool full() const {
        return true;
    }

    void sort() const {}

    igraph_real_t worstDist() const {
        return sqr_radius;
    }

    bool addPoint(igraph_real_t, igraph_int_t index) {
        if (region.contains(index)) {
            found = true;
            return false;
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

// Keeps the edges whose exclusion region contains no other point.
template<typename Region>
static igraph_error_t filter_edges(igraph_vector_int_t *edges, const igraph_matrix_t *points, igraph_real_t beta) {
    igraph_int_t available_edges = igraph_vector_int_size(edges);
    igraph_int_t added_edges = 0;
    ig_point_adaptor adaptor(points);
    igraph_int_t dim = igraph_matrix_ncol(points);
    KDTree<-1> tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();
    std::vector<igraph_real_t> centre(dim);
    Region region(points, beta);
    int iter = 0;
    for (igraph_int_t i = 0; i * 2  < available_edges; i++) {
        region.set_edge(VECTOR(*edges)[2 * i], VECTOR(*edges)[2 * i + 1]);
        bool empty = true;
        for (int k = 0; empty && k < Region::query_count; k++) {
            const igraph_real_t sqr_radius = region.query(k, centre);
            RegionSearch<Region> search(region, search_sqr_radius(centre, sqr_radius));
            tree.findNeighbors(search, centre.data());
            empty = search.size() == 0;
        }
        if (empty) {
            VECTOR(*edges)[added_edges * 2]     = VECTOR(*edges)[i * 2];
            VECTOR(*edges)[added_edges * 2 + 1] = VECTOR(*edges)[i * 2 + 1];
            added_edges += 1;
        }
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }
    IGRAPH_CHECK(igraph_vector_int_resize(edges, added_edges * 2));
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_lune_beta_skeleton
 * \brief The lune-based β-skeleton of a spatial point set.
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
 *     graph. Each row is a point; dimensionality is inferred from the column count.
 *
 * \return Error code.
 *
 * Time complexity: Around O(n^floor(d/2) log n), where n is the number of points
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
        IGRAPH_CHECK((filter_edges<LuneRegion<true>>(&potential_edges, points, beta)));
    } else {
        if (igraph_matrix_ncol(points) != 2) {
            IGRAPH_ERROR("Beta skeletons with beta < 1 are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
        }

        IGRAPH_CHECK((filter_edges<CircleIntersectionRegion>(&potential_edges, points, beta)));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_circle_beta_skeleton
 * \brief The circle-based β-skeleton of a 2D spatial point set.
 *
 * \experimental
 *
 * This function constructs the circle-based β-skeleton of a 2D spatial point set.
 *
 * </para><para>
 * A larger \p beta value results in a larger region, and a sparser graph.
 * Values of beta &lt; 1 are considerably slower.
 *
 * \param graph A pointer to the graph that will be created.
 * \param points An n-by-2 matrix containing the points that will be used to
 *    create the graph. Each row is a point.
 * \param beta A positive real value used to parameterize the graph.
 * \return Error code.
 *
 * Time complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_circle_beta_skeleton(igraph_t *graph, const igraph_matrix_t *points, igraph_real_t beta) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    if (igraph_matrix_ncol(points) != 2) {
        IGRAPH_ERROR("Circle-based beta skeletons are only supported in 2 dimensions.", IGRAPH_UNIMPLEMENTED);
    }

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, beta));
    if (beta == 1) {
        // Same test as igraph_gabriel_graph(), so that the results are identical.
        IGRAPH_CHECK(filter_edges<LuneRegion<true>>(&potential_edges, points, beta));
    } else if (beta >= 1) {
        IGRAPH_CHECK(filter_edges<CircleUnionRegion>(&potential_edges, points, beta));
    } else {
        IGRAPH_CHECK((filter_edges<CircleIntersectionRegion>(&potential_edges, points, beta)));
    }

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}

// Derived from code by Szabolcs at github.com/szhorvat/IGraphM

// Finds the smallest beta for which the closed lune of a-b contains another point.
// This is the beta at which the edge disappears from the lune-based beta skeleton.
class BetaFinder : public EdgeRegion {
    const igraph_real_t max_beta;
    igraph_real_t gabriel_band = 0;
    std::vector<igraph_real_t> centre;

    igraph_real_t smallest_beta = IGRAPH_INFINITY;
    igraph_real_t sqr_radius = 0;

    // Squared search radius around the midpoint that covers the lune for beta.
    igraph_real_t lune_search_radius(igraph_real_t beta) const {
        if (beta == IGRAPH_INFINITY) {
            return IGRAPH_INFINITY;
        }
        return search_sqr_radius(centre, 0.25 * (2 * beta - 1) * uu + boundary_band(0.25 * beta * beta * uu, true));
    }

public:
    using DistanceType = igraph_real_t;
    using IndexType = igraph_int_t;

    BetaFinder(const igraph_matrix_t *points, igraph_real_t max_beta) :
        EdgeRegion(points), max_beta(max_beta), centre(dim) { }

    void set_edge(igraph_int_t a, igraph_int_t b) {
        EdgeRegion::set_edge(a, b);
        gabriel_band = boundary_band(0.25 * uu, true);
        midpoint(centre);
        init();
    }

    void init() {
        clear();
    }

    void clear() {
        smallest_beta = IGRAPH_INFINITY;
        sqr_radius = lune_search_radius(max_beta);
    }

    bool full() const {
        return true;
    }

    size_t size() const {
        return 1;
    }

    // The beta at which point p enters the lune of a-b, see LuneRegion.
    // Returns 0 if p is in the closed Gabriel ball, i.e. a-b is not a Gabriel edge,
    // using the same test as igraph_gabriel_graph().
    igraph_real_t pointBeta(igraph_int_t p) const {
        if (p == a || p == b) {
            return IGRAPH_INFINITY;
        }
        igraph_real_t vv, ww, vw;
        local_products(p, vv, ww, vw);
        if (vw <= gabriel_band) {
            return 0;
        }
        // p is in the lune when beta >= vv / (vv - vw) and beta >= ww / (ww - vw).
        // When a denominator is not positive, p does not project strictly
        // between a and b, and is in no lune. Outside of the Gabriel ball,
        // the larger of the two bounds is the one for the endpoint closer to p.
        const igraph_real_t da = vv - vw, db = ww - vw;
        if (da <= 0 || db <= 0) {
            return IGRAPH_INFINITY;
        }
        return vv <= ww ? vv / da : ww / db;
    }

    bool addPoint(igraph_real_t, igraph_int_t index) {
        igraph_real_t beta = pointBeta(index);
        if (beta < smallest_beta && beta <= max_beta) {
            smallest_beta = beta;
            if (beta == 0) {
                return false;
            }
            sqr_radius = lune_search_radius(beta);
        }
        return true;
    }

    igraph_real_t worstDist() const {
        return sqr_radius;
    }

    void sort() {}

    const igraph_real_t *searchCentre() const {
        return centre.data();
    }

    igraph_real_t thresholdBeta() const {
        return smallest_beta;
    }
};

/**
 * \function igraph_beta_weighted_gabriel_graph
 * \brief A Gabriel graph, with edges weighted by the β value at which they disappear.
 *
 * \experimental
 *
 * This function generates a Gabriel graph, and for each edge of this graph it
 * computes the threshold β value at which the edge ceases to be part of the
 * lune-based β-skeleton. For edges that continue to be part of β-skeletons
 * for arbitrarily large β, \c IGRAPH_INFINITY is returned.
 *
 * </para><para>
 * The \p max_beta cutoff parameter controls the largest β value to consider.
 * For edges that persist above this β value, \c IGRAPH_INFINITY is returned.
 * This parameter serves to improve performance: the smaller this cutoff,
 * the faster the computation. Pass \c IGRAPH_INFINITY to use no cutoff.
 * The cutoff does not affect which edges are part of the result.
 *
 * \param graph A pointer to the graph that will be created.
 * \param weights Will contain the edge weights corresponding to the edge
 *    indices from the graph.
 * \param points A matrix containing the points that will be used.
 *    Each row is a point; dimensionality is inferred from the column count.
 *    There must be no duplicate points.
 * \param max_beta Maximum value of beta to search to, higher values will be
 *    represented as \c IGRAPH_INFINITY. Must be at least 1.
 * \return Error code.
 *
 * \sa \ref igraph_lune_beta_skeleton() or \ref igraph_circle_beta_skeleton()
 * to generate a graph with a given value of beta; \ref igraph_gabriel_graph()
 * to only generate a Gabriel graph, without edge weights.
 *
 *
 * Time complexity: Around O(n^floor(d/2) log n), where n is the number of points
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

    if (!(max_beta >= 1)) {
        IGRAPH_ERRORF("The beta cutoff must be at least 1, got %g.", IGRAPH_EINVAL, max_beta);
    }

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    igraph_int_t dim = igraph_matrix_ncol(points);
    igraph_int_t point_count = igraph_matrix_nrow(points);
    ig_point_adaptor adaptor(points);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&edges, points, 1));
    igraph_int_t edge_count = igraph_vector_int_size(&edges) / 2;

    IGRAPH_CHECK(igraph_vector_resize(weights, edge_count));

    KDTree<-1> tree(dim, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

    BetaFinder finder(points, max_beta);
    int iter = 0;
    for (igraph_int_t i = 0; i < edge_count; i++) {
        finder.set_edge(VECTOR(edges)[2 * i], VECTOR(edges)[2 * i + 1]);
        tree.findNeighbors(finder, finder.searchCentre());
        VECTOR(*weights)[i] = finder.thresholdBeta();
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
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
    IGRAPH_FINALLY_CLEAN(1);

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
 * computing the Gabriel graph of point sets of arbitrary dimension.
 *
 * </para><para>
 * The Gabriel graph is a special case of lune-based and circle-based β-skeletons
 * with <code>β=1</code>.
 *
 * \param graph A pointer to the graph to be created.
 * \param points The point set that will be used. Each row is a point,
 *    dimensionality is inferred from column count.
 *
 * \return Error code.
 * \sa The Gabriel graph is a special case of
 *     \ref igraph_lune_beta_skeleton() and \ref igraph_circle_beta_skeleton()
 *     where <code>β = 1</code>.
 *
 * Time complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_gabriel_graph(igraph_t *graph, const igraph_matrix_t *points) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, 1));

    IGRAPH_CHECK((filter_edges<LuneRegion<true>>(&potential_edges, points, 1)));

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
 * \sa \ref igraph_lune_beta_skeleton() to compute the lune-based β-skeleton
 * for <code>β = 2</code> or other β values.
 *
 * Time complexity: Around O(n^floor(d/2) log n), where n is the number of points
 * and d is the dimensionality of the point set.
 *
 */
igraph_error_t igraph_relative_neighborhood_graph(igraph_t *graph, const igraph_matrix_t *points) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_vector_int_t potential_edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&potential_edges, 0);

    IGRAPH_CHECK(beta_skeleton_edge_superset(&potential_edges, points, 1));

    IGRAPH_CHECK((filter_edges<LuneRegion<false>>(&potential_edges, points, 2)));

    IGRAPH_CHECK(igraph_create(graph, &potential_edges, igraph_matrix_nrow(points), false));

    igraph_vector_int_destroy(&potential_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}
