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
#include "igraph_interface.h"
#include "igraph_spatial.h"

#include "igraph_conversion.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "core/exceptions.h"

#include "nanoflann/nanoflann.hpp"

#include <vector>

class ig_point_adaptor {
    const igraph_matrix_t *points;
    const igraph_integer_t point_count;

public:
    explicit ig_point_adaptor(const igraph_matrix_t *points) :
        points(points), point_count(igraph_matrix_nrow(points)) { }

    size_t kdtree_get_point_count() const {
        return point_count;
    }

    igraph_real_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        return MATRIX(*points, idx, dim);
    }

    template <typename BoundingBox>
    bool kdtree_get_bbox(BoundingBox &bb) const {
        IGRAPH_UNUSED(bb);
        return false; // indicates that it should use default
    }
};


class GraphBuildingResultSet {
    igraph_integer_t current_vertex = 0;
    size_t added_count = 0;
    const igraph_real_t max_distance;
    const igraph_integer_t max_neighbors;

public:
    std::vector<igraph_integer_t> neighbors;
    std::vector<igraph_real_t> distances;

    using DistanceType = igraph_real_t;

    GraphBuildingResultSet(const igraph_integer_t max_neighbors, const igraph_real_t max_distance) :
        max_distance(max_distance),
        max_neighbors(max_neighbors),
        neighbors(max_neighbors),
        distances(max_neighbors) { }

    bool addPoint(const igraph_real_t distance, const igraph_integer_t index) {
        size_t i;

        if (index == current_vertex) {
            return true;
        }

        for (i = added_count; i > 0; i--) {
            if ((distances[i-1] > distance /*|| dists[i-1] == distance && index < edges[i-1]*/ )) {
                if (i < max_neighbors) {
                    distances[i] = distances[i-1];
                    neighbors[i] = neighbors[i-1];
                }
            } else {
                break;
            }
        }
        if (i < max_neighbors) {
            neighbors[i] = index;
            distances[i] = distance;
        }
        if (added_count != max_neighbors) {
            added_count++;
        }
        return true;
    }

    void reset(igraph_integer_t current_vertex_) {
        added_count = 0;
        current_vertex = current_vertex_;
    }

    void sort() { }

    size_t size() const {
        return added_count;
    }

    bool full() const {
        return added_count == max_neighbors;
    }

    bool empty() const {
        return added_count == 0;
    }

    igraph_real_t worstDist() const {
        if (added_count < max_neighbors || added_count == 0) {
            return  max_distance;
        }
        return distances[added_count-1];
    }
};

template <typename Metric, igraph_integer_t Dimension>
static igraph_error_t neighbor_helper(
        igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_integer_t neighbors,
        igraph_real_t cutoff,
        igraph_integer_t dimension,
        igraph_bool_t directed) {

    const igraph_integer_t point_count = igraph_matrix_nrow(points);

    using kdTree = nanoflann::KDTreeSingleIndexAdaptor<Metric, ig_point_adaptor, Dimension>;

    ig_point_adaptor adaptor(points);

    kdTree tree(dimension, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));

    tree.buildIndex();

    igraph_vector_t current_point;
    IGRAPH_VECTOR_INIT_FINALLY(&current_point, dimension);

    igraph_integer_t neighbor_count = neighbors >= 0 ? neighbors : point_count;

    GraphBuildingResultSet results(neighbor_count, cutoff);
    std::vector<igraph_integer_t> edges;
    for (igraph_integer_t i = 0; i < point_count; i++) {

        results.reset(i);
        IGRAPH_CHECK(igraph_matrix_get_row(points, &current_point, i));

        tree.findNeighbors(results, VECTOR(current_point), nanoflann::SearchParameters(0, false));
        for (igraph_integer_t j = 0; j < results.size(); j++) {
            edges.push_back(i);
            edges.push_back(results.neighbors[j]);
        }
    }

    igraph_vector_destroy(&current_point);
    IGRAPH_FINALLY_CLEAN(1);

    // Overflow check, ensures that edges.size() is not too large for igraph vectors.
    if (edges.size() > IGRAPH_INTEGER_MAX) {
        IGRAPH_ERROR("Too many edges.", IGRAPH_EOVERFLOW);
    }

    igraph_vector_int_t edge_view;
    igraph_vector_int_view(&edge_view, edges.data(), edges.size());
    IGRAPH_CHECK(igraph_create(graph, &edge_view, point_count, true));

    if (! directed) {
        IGRAPH_CHECK(igraph_to_undirected(graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL));
    }

    return IGRAPH_SUCCESS;
}


template <typename Metric>
static igraph_error_t dimension_dispatcher(
        igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_integer_t neighbors,
        igraph_real_t cutoff,
        igraph_integer_t dimension,
        igraph_bool_t directed) {

    switch (dimension) {
    case 0:
        IGRAPH_ERROR("0-dimensional points are not supported.", IGRAPH_EINVAL);
    case 1:
        return neighbor_helper<Metric, 1>(graph, points, neighbors, cutoff, dimension, directed);
    case 2:
        return neighbor_helper<Metric, 2>(graph, points, neighbors, cutoff, dimension, directed);
    case 3:
        return neighbor_helper<Metric, 3>(graph, points, neighbors, cutoff, dimension, directed);
    default:
        return neighbor_helper<Metric, -1>(graph, points, neighbors, cutoff, dimension, directed);
    }
}


/**
 * \function igraph_nearest_neighbor_graph
 * \brief Computes the nearest neighbor graph for a spatial point set.
 *
 * This function constructs the \p k nearest neighbor graph of a given point
 * set. Each point is connected to at most \p k spatial neighbors within a
 * radius of \p cutoff.
 *
 * \param graph A pointer to the graph that will be created.
 * \param points A matrix containing the points that will be used to create
 *    the graph. Each row is a point, dimensionality is inferred from the
 *    column count.
 * \param metric The distance metric to use. See \ref igraph_metric_t.
 * \param neighbors How many neighbors will be added for each vertex, set to
 *    a negative value to ignore.
 * \param cutoff Maximum distance at which connections will be made, set to a
 *    negative value or \c IGRAPH_INFINITY to ignore.
 * \param directed Whether to create a directed graph.
 * \return Error code.
 */
igraph_error_t igraph_nearest_neighbor_graph(igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_metric_t metric,
        igraph_integer_t k,
        igraph_real_t cutoff,
        igraph_bool_t directed) {

    const igraph_integer_t dimension = igraph_matrix_ncol(points);

    // Negative cutoff values signify that no cutoff should be used.
    cutoff = cutoff >= 0 ? cutoff : IGRAPH_INFINITY;

    // Handle null graph separately.
    // The number of matrix columns is not meaningful when there are zero rows.
    // This prevents throwing an error when both the column and the row counts are zero.
    if (igraph_matrix_nrow(points) == 0) {
        return igraph_empty(graph, 0, directed);
    }

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    switch (metric) {
    case IGRAPH_METRIC_L2:
        return dimension_dispatcher<nanoflann::L2_Adaptor<igraph_real_t, ig_point_adaptor> > (
                   graph,
                   points,
                   k,
                   cutoff * cutoff, // L2 uses square distances, so adjust for that here.
                   dimension,
                   directed);
    case IGRAPH_METRIC_L1:
            return dimension_dispatcher<nanoflann::L1_Adaptor<igraph_real_t, ig_point_adaptor> > (
                    graph,
                    points,
                    k,
                    cutoff,
                    dimension,
                    directed);
    default:
        IGRAPH_ERROR("Invalid metric.", IGRAPH_EINVAL);
    }

    IGRAPH_HANDLE_EXCEPTIONS_END;
}
