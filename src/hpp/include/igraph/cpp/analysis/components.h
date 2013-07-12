/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ANALYSIS_COMPONENTS_H
#define IGRAPHPP_ANALYSIS_COMPONENTS_H

#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>

namespace igraph {

class Graph;

/// Finds the (weakly or strongly) connected components of a graph
void clusters(const Graph& graph, Vector* membership,
        Vector* csize = 0, integer_t* no = 0,
        Connectedness connectedness = IGRAPH_WEAK);

}         // end of namespace

#endif    // IGRAPHPP_ANALYSIS_COMPONENTS_H

