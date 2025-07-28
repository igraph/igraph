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

public:
    explicit ig_point_adaptor(const igraph_matrix_t *points) :
        points(points) { }

    size_t kdtree_get_point_count() const {
        return igraph_matrix_nrow(points);
    }

    igraph_real_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        return MATRIX(*points, idx, dim);
    }

    // indicates that it should use default
    template <typename BoundingBox>
    bool kdtree_get_bbox(BoundingBox &bb) const {
        IGRAPH_UNUSED(bb);
        return false;
    }
};


class GraphBuildingResultSet {
    igraph_integer_t current_vertex = 0;
    igraph_integer_t current_added = 0;
    igraph_real_t max_distance;
    igraph_integer_t max_neighbors;
    igraph_integer_t *edges;
    igraph_real_t *dists;

public:
    using DistanceType = igraph_real_t;

    GraphBuildingResultSet(const igraph_integer_t neighbors, const igraph_real_t distance) :
        max_distance(distance),
        max_neighbors(neighbors) {}

    bool addPoint(const igraph_real_t distance, const igraph_integer_t index) {
        igraph_integer_t i;

        if (index == current_vertex) {
            return true;
        }

        for (i = current_added; i > 0; i--) {
            if ((dists[i-1] > distance /*|| dists[i-1] == distance && index < edges[i-1]*/ )) {
                    if (i < max_neighbors) {
                    dists[i] = dists[i-1];
                    edges[i] = edges[i-1];
                }
            } else {
                break;
            }
        }
        if (i < max_neighbors) {
            edges[i] = index;
            dists[i] = distance;
        }
        if (current_added != max_neighbors) {
            current_added++;
        }
        return true;
    }

    void init(igraph_real_t *dists_, igraph_integer_t *edges_, igraph_integer_t current_vertex_) {
        edges = edges_;
        dists = dists_;
        current_added = 0;
        current_vertex = current_vertex_;
    }

    void sort() { }

    igraph_integer_t size() const {
        return current_added;
    }

    bool full() const {
        return current_added == max_neighbors;
    }

    bool empty() const {
        return current_added == 0;
    }

    igraph_real_t worstDist() const {
        if (current_added < max_neighbors || current_added == 0) {
            return  max_distance;
        }
        return dists[current_added-1];
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

    kdTree tree(dimension, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(1));

    tree.buildIndex();

    igraph_vector_t current_point;
    IGRAPH_VECTOR_INIT_FINALLY(&current_point, dimension);

    igraph_integer_t neighbor_count = neighbors >= 0 ? neighbors : point_count;
    cutoff = cutoff >= 0 ? cutoff : INFINITY;

    using resultClass = GraphBuildingResultSet;
    resultClass results(neighbor_count, cutoff);
    std::vector<igraph_integer_t> neighbor_set(neighbor_count);
    std::vector<igraph_real_t>    distances(neighbor_count);
    std::vector<igraph_integer_t> edges;
    for (igraph_integer_t i = 0; i < point_count; i++) {

        std::fill(neighbor_set.begin(), neighbor_set.end(), 0);
        std::fill(distances.begin(), distances.end(), 0);

        results.init(distances.data(), neighbor_set.data(), i);
        IGRAPH_CHECK(igraph_matrix_get_row(points, &current_point, i));

        tree.findNeighbors(results, VECTOR(current_point), nanoflann::SearchParameters(0, false));
        for (igraph_integer_t j = 0; j < results.size(); j++) {
            edges.push_back(i);
            edges.push_back(neighbor_set[j]);
        }
    }

    igraph_vector_destroy(&current_point);
    IGRAPH_FINALLY_CLEAN(1);

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
        IGRAPH_ERROR("0-dimensional points are not supported", IGRAPH_EINVAL);
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
 *    negative value or \c INFINITY to ignore.
 * \param directed Whether to create a directed graph.
 * \return Error code.
 */
igraph_error_t igraph_nearest_neighbor_graph(igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_metric_t metric,
        igraph_integer_t k,
        igraph_real_t cutoff,
        igraph_bool_t directed) {

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    const igraph_integer_t dimension = igraph_matrix_ncol(points);
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
