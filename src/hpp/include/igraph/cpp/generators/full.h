/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_FULL_H
#define IGRAPHPP_GENERATORS_FULL_H

#include <memory>

namespace igraph {

class Graph;

/// Generates a full graph with the given number of nodes
std::auto_ptr<Graph> full(integer_t nodes, bool directed = false,
        bool loops = false);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_FULL_H
