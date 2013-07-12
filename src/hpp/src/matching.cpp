/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/graph.h>
#include <igraph/cpp/vector.h>
#include <igraph/cpp/vector_bool.h>
#include <igraph/cpp/vector_long.h>
#include <igraph/igraph_matching.h>

namespace igraph {

void maximum_bipartite_matching(const Graph& graph, const VectorBool& types,
        integer_t* matching_size, real_t* matching_weight, VectorLong* matching,
        const Vector* weights, real_t eps) {
    IGRAPH_TRY(igraph_maximum_bipartite_matching(graph.c_graph(), types.c_vector(),
                matching_size, matching_weight,
                matching ? matching->c_vector() : 0,
                weights ? weights->c_vector() : 0,
                eps
    ));
}

}          // end of namespace

