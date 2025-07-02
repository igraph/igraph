
#include "cliques/cliquer/graph.h"
#include "igraph_operators.h"
#include "igraph_spatial.h"
#include "igraph_datatype.h"
#include "igraph_matrix.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "nanoflann/nanoflann.hpp"
#include <cstddef>
#include <cstdint>
#include <string.h>

class igraph_point_adaptor {
    igraph_integer_t dimension;
    const igraph_matrix_t *points;
    size_t kdtree_get_point_count() const {
        return igraph_matrix_nrow(points);
    }
    igraph_real_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        return igraph_matrix_get(points, idx, dim);
    }

    // indicates that it should use default
    bool kdtree_get_bbox(void* bb) {
        return false;
    }

public:
    igraph_point_adaptor(const igraph_matrix_t *points) {
        this -> dimension = igraph_matrix_ncol(points);
        this -> points = points;
    }
};

class L2_igraph_adaptor {
public:
    using ElementType = igraph_real_t;
    using DistanceType = nanoflann::metric_L2;

} ;

class GraphBuildingResultSet {
    private:
        igraph_integer_t current_vertex = 0;
        igraph_t *graph;
    public:
        using DistanceType = igraph_real_t;
        GraphBuildingResultSet(igraph_t * graph, igraph_integer_t neighbors, igraph_real_t distance) {
           this->graph = graph;
        }
        bool addPoint(igraph_real_t distance, igraph_integer_t index) {
            igraph_add_edge(graph, current_vertex, index);
            return true;
        }
        void select_vertex(igraph_integer_t i) {
            current_vertex = i;
        }
};

extern "C" igraph_error_t igraph_nearest_neighbor_graph(igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_metric_t metric,
        igraph_integer_t neighbors,
        igraph_real_t cutoff) {

    igraph_integer_t dimensionality = igraph_matrix_ncol(points);

    igraph_point_adaptor adaptor(points);

    igraph_integer_t point_count = igraph_matrix_nrow(points);

    using my_tree = nanoflann::KDTreeSingleIndexAdaptor<L2_igraph_adaptor, igraph_point_adaptor, 0, uint32_t>;

    my_tree tree(dimensionality, adaptor, 10);

    tree.buildIndex();

    igraph_vector_t current_point;
    IGRAPH_VECTOR_INIT_FINALLY(&current_point, 0);

    GraphBuildingResultSet results(graph, 2, 100);


    for (igraph_integer_t i = 0; i < point_count; i++) {
        IGRAPH_CHECK(igraph_matrix_get_row(points, &current_point, i));
        results.select_vertex(i);
        tree.findNeighbors(results, VECTOR(current_point), nanoflann::SearchParameters(0, false));
    }

    igraph_vector_destroy(&current_point);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
