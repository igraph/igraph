#include "igraph_constructors.h"
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

    kdTree tree(dimension, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(1));

    tree.buildIndex();

    igraph_vector_t current_point;
    IGRAPH_VECTOR_INIT_FINALLY(&current_point, 0);

    igraph_integer_t neighbor_count = neighbors > 0 ? neighbors : point_count;

    using resultClass = nanoflann::RKNNResultSet<igraph_real_t, igraph_integer_t, igraph_integer_t>;
    resultClass results(neighbor_count, cutoff);
    std::vector<igraph_integer_t> neighbor_set(neighbor_count);
    std::vector<igraph_real_t>    distances(neighbor_count);
    std::vector<igraph_integer_t> edges;
    for (igraph_integer_t i = 0; i < point_count; i++) {

        std::fill(neighbor_set.begin(), neighbor_set.end(), 0);
        std::fill(distances.begin(), distances.end(), 0);

        results.init(neighbor_set.data(), distances.data());
        IGRAPH_CHECK(igraph_matrix_get_row(points, &current_point, i));

        tree.findNeighbors(results, VECTOR(current_point), nanoflann::SearchParameters(0, false));
        printf("beep %li, %li\n", results.size(), neighbor_set.size());
        for (igraph_integer_t j = 0; j < results.size(); j++) {
            printf("adding %li->%li\n", i, j);
            edges.push_back(i);
            edges.push_back(j);
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
