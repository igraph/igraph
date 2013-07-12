/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_IO_LGL_H
#define IGRAPHPP_IO_LGL_H

#include <igraph/cpp/graph.h>

namespace igraph {

/// Reads a graph from an LGL file
Graph read_lgl(FILE* instream, bool names=true,
        AddWeights weights = IGRAPH_ADD_WEIGHTS_IF_PRESENT,
        bool directed=true);

/// Writes a graph to an LGL file
void write_lgl(const Graph& graph, FILE* outstream, const std::string& names="name",
        const std::string& weights="weight", bool isolates=true);

}         // end of namespace

#endif    // IGRAPHPP_IO_LGL_H

