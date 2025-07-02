
#include "igraph_spatial.h"
#include "igraph_datatype.h"
#include "igraph_matrix.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "nanoflann/nanoflann.hpp"
#include <cstdint>

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

extern "C" igraph_error_t igraph_nearest_neighbor_graph(igraph_t *graph,
        const igraph_matrix_t *points,
        igraph_metric_t metric,
        igraph_integer_t neighbors,
        igraph_real_t cutoff) {

    igraph_integer_t dimensionality = igraph_matrix_ncol(points);

    igraph_point_adaptor adaptor(points);

    using my_tree = nanoflann::KDTreeSingleIndexAdaptor<L2_igraph_adaptor, igraph_point_adaptor, 0, uint32_t>;

    my_tree tree(dimensionality, adaptor, 10);


    tree.buildIndex();

    return IGRAPH_SUCCESS;
}
