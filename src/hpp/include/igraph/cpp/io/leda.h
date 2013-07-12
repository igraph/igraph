/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_IO_LEDA_H
#define IGRAPHPP_IO_LEDA_H

#include <igraph/cpp/graph.h>

namespace igraph {

/// Writes the graph in LEDA format to the given file
void write_leda(const Graph& graph, FILE* outstream,
        const std::string& vertex_attr_name = "",
        const std::string& edge_attr_name = "");

}         // end of namespace

#endif    // IGRAPHPP_IO_LEDA_H

