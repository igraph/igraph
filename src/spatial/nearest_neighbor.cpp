/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
u
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_constructors.h"
#include "igraph_spatial.h"

#include "igraph_interface.h"
#include "igraph_matrix.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"

#include "nanoflann/nanoflann.hpp"
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string.h>
#include <vector>

class igraph_point_adaptor {
    igraph_integer_t dimension;
    const igraph_matrix_t *points;

public:
    size_t kdtree_get_point_count() const {
        return igraph_matrix_nrow(points);
    }
    igraph_real_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        return igraph_matrix_get(points, idx, dim);
    }

    // indicates that it should use default
    template <typename BoundingBox>
    bool kdtree_get_bbox(BoundingBox &bb) const {
        return false;
    }

public:
    igraph_point_adaptor(const igraph_matrix_t *points) {
        this -> dimension = igraph_matrix_ncol(points);
        this -> points = points;
    }
};
/*
class L2_igraph_adaptor {
public:
    using ElementType = igraph_real_t;
    using DistanceType = nanoflann::metric_L2;
    L2_igraph_adaptor(const igraph_point_adaptor& point) {

    }
} ;
*/

class GraphBuildingResultSet {
private:
    igraph_integer_t current_vertex = 0;
    igraph_integer_t current_added = 0;
    igraph_real_t max_distance;
    igraph_integer_t max_neighbors;
    igraph_integer_t *edges;
    igraph_real_t *dists;
public:
    using DistanceType = igraph_real_t;
    GraphBuildingResultSet(igraph_integer_t neighbors, igraph_real_t distance) :
        max_distance(distance),
        max_neighbors(neighbors) {}

    bool addPoint(igraph_real_t distance, igraph_integer_t index) {
        if (index == current_vertex) {
            return true;
        }
        igraph_integer_t i;
        for (i = current_added; i > 0; i--) {
            if ((dists[i-1] < distance || dists[i-1] == distance && index < edges[i-1] ) && i < max_neighbors) {
                dists[i] = dists[i-1];
                edges[i] = edges[i-1];
            } else { break; }
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

        void init(igraph_real_t *_dists, igraph_integer_t* _edges, igraph_integer_t current_vertex) {
            edges = _edges;
            dists = _dists;
            current_added = 0;
            this-> current_vertex = current_vertex;
        }

        void sort () {}

        igraph_integer_t size() {
            return current_added;
        }

    bool full () const {
        return current_added == max_neighbors;
    }
        bool empty() const {
            return current_added == 0;
        }

    igraph_real_t worstDist() {
        if (current_added < max_neighbors || current_added == 0) {
            return  max_distance;
        }
        return dists[current_added-1];
    }

};

template <typename Metric, uint32_t Dimension>
static igraph_error_t neighbor_helper(
    igraph_t *graph,
    const igraph_matrix_t *points,
    igraph_integer_t neighbors,
    igraph_real_t cutoff,
    igraph_integer_t dimension) {

    igraph_integer_t point_count = igraph_matrix_nrow(points);

    using kdTree = nanoflann::KDTreeSingleIndexAdaptor<Metric, igraph_point_adaptor, Dimension>;

    igraph_point_adaptor adaptor(points);

    kdTree tree(dimension, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(4));

    tree.buildIndex();

    igraph_vector_t current_point;
    IGRAPH_VECTOR_INIT_FINALLY(&current_point, 0);

    igraph_integer_t neighbor_count = neighbors > 0 ? neighbors : point_count;

    //using resultClass = nanoflann::RKNNResultSet<igraph_real_t, igraph_integer_t, igraph_integer_t>;
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
    igraph_vector_int_t edge_view;
    igraph_vector_int_view(&edge_view, edges.data(), edges.size());
    igraph_create(graph, &edge_view, point_count, true);

    //igraph_vector_int_destroy(&edge_view);

    igraph_vector_destroy(&current_point);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}


template <typename Metric>
static igraph_error_t dimension_dispatcher(
    igraph_t *graph,
    const igraph_matrix_t *points,
    igraph_integer_t neighbors,
    igraph_real_t cutoff,
    igraph_integer_t dimension) {
    switch (dimension) {
    case 2:  return neighbor_helper<Metric, 2>(graph, points, neighbors, cutoff, dimension);
    case 3:  return neighbor_helper<Metric, 3>(graph, points, neighbors, cutoff, dimension);
    default: return neighbor_helper<Metric, 0>(graph, points, neighbors, cutoff, dimension);
    }
}

igraph_error_t igraph_nearest_neighbor_graph(igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_metric_t metric,
        igraph_integer_t neighbors,
        igraph_real_t cutoff) {
    igraph_integer_t dimension = igraph_matrix_ncol(points);
    switch (metric) {
    case IGRAPH_METRIC_L2 :
        return dimension_dispatcher<nanoflann::L2_Adaptor<igraph_real_t, igraph_point_adaptor> > (
                   graph,
                   points,
                   neighbors,
                   cutoff * cutoff, // L2 uses square distances, so adjust for that here.
                   dimension);
    default : IGRAPH_ERROR("Metic type not implemented.", IGRAPH_UNIMPLEMENTED);
    }
}
