/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/igraph_structural.h>
#include <igraph/cpp/edge_selector.h>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/vector.h>
#include <igraph/cpp/vector_bool.h>

namespace igraph {

Vector count_multiple(const EdgeSelector& es) {
    Vector result;
    IGRAPH_TRY(igraph_count_multiple(es.getGraph()->c_graph(),
                result.c_vector(), *es.c_es()));
    return result;
}

bool has_multiple(const Graph& graph) {
    igraph_bool_t result;
    IGRAPH_TRY(igraph_has_multiple(graph.c_graph(), &result));
    return result;
}

VectorBool is_multiple(const EdgeSelector& es) {
    VectorBool result;
    IGRAPH_TRY(igraph_is_multiple(es.getGraph()->c_graph(),
                result.c_vector(), *es.c_es()));
    return result;
}

bool is_simple(const Graph& graph) {
    igraph_bool_t result;
    IGRAPH_TRY(igraph_is_simple(graph.c_graph(), &result));
    return result;
}

}         // end of namespaces
