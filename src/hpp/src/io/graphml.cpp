/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/error.h>
#include <igraph/cpp/io/graphml.h>
#include <memory>

namespace igraph {

Graph read_graphml(FILE* instream, int index) {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_read_graph_graphml(result.get(), instream, index));
    return Graph(result.release());
}

void write_graphml(const Graph& graph, FILE* outstream) {
    IGRAPH_TRY(igraph_write_graph_graphml(graph.c_graph(), outstream));
}

}         // end of namespace

