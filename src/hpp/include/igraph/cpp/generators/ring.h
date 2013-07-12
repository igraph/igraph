/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_RING_H
#define IGRAPHPP_GENERATORS_RING_H

#include <memory>

namespace igraph {

class Graph;

/// Generates a ring graph
std::auto_ptr<Graph> ring(integer_t n, bool directed = false, bool mutual = false,
        bool circular = true);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_RING_H
