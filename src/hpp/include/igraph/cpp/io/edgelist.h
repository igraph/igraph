/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_IO_EDGELIST_H
#define IGRAPHPP_IO_EDGELIST_H

#include <igraph/cpp/graph.h>

namespace igraph {

/// Reads a graph from an edge list file
Graph read_edgelist(FILE* instream, integer_t n=0, bool directed=true);

/// Writes the edge list of the graph to the given file
void write_edgelist(const Graph& graph, FILE* outstream);

}         // end of namespace

#endif    // IGRAPHPP_IO_EDGELIST_H

