#include "igraph_interface.h"
#include "igraph_operators.h"
#include "igraph_spatial.h"
#include "igraph_matrix.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "nanoflann/nanoflann.hpp"
#include <cstddef>
#include <cstdint>
#include <string.h>

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
    igraph_integer_t max_neighbors;
    igraph_real_t max_distance;
    igraph_vector_int_t edges;
    igraph_t *graph;
public:
    using DistanceType = igraph_real_t;
    GraphBuildingResultSet(igraph_t * graph, igraph_integer_t neighbors, igraph_real_t distance) : max_distance(distance), max_neighbors(neighbors) {
        this->graph = graph;
    }
    igraph_error_t initialize() {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        return IGRAPH_SUCCESS;
    }
    bool addPoint(igraph_real_t distance, igraph_integer_t index) {
        if (index == current_vertex) {
            return true;
        }
        igraph_vector_int_push_back(&edges, current_vertex);
        igraph_vector_int_push_back(&edges, index);
        if (++current_added == max_neighbors) {
            return false;
        }
        return true;
    }
    ~GraphBuildingResultSet() {
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }
    void select_vertex(igraph_integer_t i) {
        current_vertex = i;
        current_added = 0;
    }
    void sort () {

    }
    bool full () {
        return false;
    }
    igraph_real_t worstDist() {
        return max_distance;
    }
    igraph_error_t dumpGraph() {
        IGRAPH_CHECK(igraph_add_edges(graph, &edges, NULL));
        return IGRAPH_SUCCESS;
    }
};

igraph_error_t igraph_nearest_neighbor_graph(igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_metric_t metric,
        igraph_integer_t neighbors,
        igraph_real_t cutoff) {

    igraph_integer_t dimensionality = igraph_matrix_ncol(points);

    igraph_integer_t point_count = igraph_matrix_nrow(points);

    using kdTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Adaptor<igraph_real_t, igraph_point_adaptor>, igraph_point_adaptor>;

    igraph_empty(graph, point_count, false);

    igraph_point_adaptor adaptor(points);

    kdTree tree(dimensionality, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(1));

    tree.buildIndex();

    igraph_vector_t current_point;
    IGRAPH_VECTOR_INIT_FINALLY(&current_point, 0);

    GraphBuildingResultSet results(graph, neighbors, cutoff);
    IGRAPH_CHECK(results.initialize());

    for (igraph_integer_t i = 0; i < point_count; i++) {
        IGRAPH_CHECK(igraph_matrix_get_row(points, &current_point, i));
        results.select_vertex(i);
        tree.findNeighbors(results, VECTOR(current_point), nanoflann::SearchParameters(0, false));
    }

    IGRAPH_CHECK(results.dumpGraph());

    igraph_vector_destroy(&current_point);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
