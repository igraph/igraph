/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/io/leda.h>

namespace igraph {

void write_leda(const Graph& graph, FILE* outstream,
        const std::string& vertex_attr_name,
        const std::string& edge_attr_name) {
    const char *vattr, *eattr;

    vattr = vertex_attr_name.length() ? vertex_attr_name.c_str() : 0;
    eattr = edge_attr_name.length() ? edge_attr_name.c_str() : 0;

    IGRAPH_TRY(igraph_write_graph_leda(graph.c_graph(), outstream, vattr, eattr));
}

}         // end of namespace

