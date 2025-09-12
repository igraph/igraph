#ifndef SPATIAL_NANOFLANN_INTERNAL_H
#define SPATIAL_NANOFLANN_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_types.h"

#include "igraph_matrix.h"

#include "nanoflann/nanoflann.hpp"

#include <vector>


class ig_point_adaptor {
    const igraph_matrix_t *points;
    const igraph_int_t point_count;

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
    igraph_int_t added_count = 0;
    const igraph_real_t max_distance;
    const igraph_int_t max_neighbors;

public:
    igraph_int_t current_vertex = 0;
    std::vector<igraph_int_t> neighbors;
    std::vector<igraph_real_t> distances;

    using DistanceType = igraph_real_t;
    using IndexType = igraph_int_t;

    GraphBuildingResultSet(const igraph_int_t max_neighbors, const igraph_real_t max_distance) :
        max_distance(max_distance),
        max_neighbors(max_neighbors),
        neighbors(0),
        distances(0) { }

    bool addPoint(const igraph_real_t distance, const igraph_int_t index) {
        igraph_int_t i;

        if (index == current_vertex) {
            return true;
        }

        for (i = added_count; i > 0; i--) {
            // TODO: Stabilize result in case of multiple points at exactly the same distance?
            // See NANOFLANN_FIRST_MATCH in RKNNResultSet in nanoflann.hpp for reference.
            if (distances[i - 1] > distance) {
                registerPoint(i, neighbors[i - 1], distances[i - 1]);
            } else {
                break;
            }
        }
        if (i < max_neighbors) {
            registerPoint(i, index, distance);
        }
        if (added_count != max_neighbors) {
            added_count++;
        }
        return true;
    }

    void registerPoint(igraph_int_t where, igraph_int_t index, igraph_real_t distance) {
        if (neighbors.size() == where) {
            neighbors.push_back(index);
            distances.push_back(distance);
        } else  {
            neighbors[where] = index;
            distances[where] = distance;
        }

    }

    void reset(igraph_int_t current_vertex_) {
        added_count = 0;
        current_vertex = current_vertex_;
    }

    // Never called. Necessary to conform to the interface.
    void sort() { }

    igraph_int_t size() const {
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
        return distances[added_count - 1];
    }
};


#endif // SPATIAL_NANOFLANN_INTERNAL_H
