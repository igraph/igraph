/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/error.h>
#include <igraph/cpp/io/edgelist.h>
#include <memory>

namespace igraph {

Graph read_edgelist(FILE* instream, integer_t n, bool directed) {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_read_graph_edgelist(result.get(), instream, n, directed));
    return Graph(result.release());
}

void write_edgelist(const Graph& graph, FILE* outstream) {
    IGRAPH_TRY(igraph_write_graph_edgelist(graph.c_graph(), outstream));
}

}         // end of namespace

