/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/graph.h>
#include <igraph/cpp/io/lgl.h>
#include <memory>

namespace igraph {

Graph read_lgl(FILE* instream, bool names, AddWeights weights, bool directed) {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_read_graph_lgl(result.get(), instream, names, weights, directed));
    return Graph(result.release());
}

void write_lgl(const Graph& graph, FILE* outstream, const std::string& names,
        const std::string& weights, bool isolates) {
    const char* names_str = names.length() > 0 ? names.c_str() : 0;
    const char* weights_str = weights.length() > 0 ? weights.c_str() : 0;

    IGRAPH_TRY(igraph_write_graph_lgl(graph.c_graph(), outstream,
                names_str, weights_str, isolates));
}

}         // end of namespace

