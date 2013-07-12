/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_LINE_GRAPH_H
#define IGRAPHPP_GENERATORS_LINE_GRAPH_H

namespace igraph {

class Graph;

/// Returns the line graph of the given graph
std::auto_ptr<Graph> line_graph(const Graph& graph);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_LINE_GRAPH_H
