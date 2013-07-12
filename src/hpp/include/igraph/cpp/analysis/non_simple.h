/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ANALYSIS_NON_SIMPLE_H
#define IGRAPHPP_ANALYSIS_NON_SIMPLE_H

#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>

namespace igraph {

class EdgeSelector;
class Graph;

/// For each edge in the edge selector, tells its multiplicity
Vector count_multiple(const EdgeSelector& es);

/// Decides whether the input graph has multiple edges
bool has_multiple(const Graph& graph);

/// For each edge in the edge selector, tells whether it is multiple
VectorBool is_multiple(const EdgeSelector& es);

/// Decides whether the input graph is a simple graph
bool is_simple(const Graph& graph);

}         // end of namespace

#endif    // IGRAPHPP_ANALYSIS_NON_SIMPLE_H

