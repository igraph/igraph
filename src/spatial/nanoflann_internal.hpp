#ifndef NANOFLANN_INTERNAL_H_
#define NANOFLANN_INTERNAL_H_

#include "igraph_decls.h"
#include "igraph_types.h"

#include "igraph_matrix.h"

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
    igraph_integer_t added_count = 0;
    const igraph_real_t max_distance;
    const igraph_integer_t max_neighbors;

public:
    igraph_integer_t current_vertex = 0;
    std::vector<igraph_integer_t> neighbors;
    std::vector<igraph_real_t> distances;

    using DistanceType = igraph_real_t;
    using IndexType = igraph_integer_t;

    GraphBuildingResultSet(const igraph_integer_t max_neighbors, const igraph_real_t max_distance) :
        max_distance(max_distance),
        max_neighbors(max_neighbors),
        neighbors(1),
        distances(1) { }

    bool addPoint(const igraph_real_t distance, const igraph_integer_t index) {
        igraph_integer_t i;

        if (index == current_vertex) {
            return true;
        }

        for (i = added_count; i > 0; i--) {
            // TODO: Stabilize result in case of multiple points at exactly the same distance?
            // See NANOFLANN_FIRST_MATCH in RKNNResultSet in nanoflann.hpp for reference.
            if (distances[i - 1] > distance) {
                if (i < max_neighbors) {
                    distances[i] = distances[i - 1];
                    neighbors[i] = neighbors[i - 1];
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
            if (added_count >= neighbors.size()) {
                neighbors.push_back(0); // always keep space to add new point.
                distances.push_back(0);
            }
        }
        return true;
    }

    void reset(igraph_integer_t current_vertex_) {
        added_count = 0;
        current_vertex = current_vertex_;
    }

    // Never called. Necessary to conform to the interface.
    void sort() { }

    igraph_integer_t size() const {
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


#endif // NANOFLANN_INTERNAL_H_
