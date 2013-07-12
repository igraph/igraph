/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_GRG_H
#define IGRAPHPP_GENERATORS_GRG_H

#include <memory>

namespace igraph {

class Graph;

/// Generates a geometric random graph
std::auto_ptr<Graph> grg_game(integer_t nodes, real_t radius, bool torus = false,
        Vector* x = 0, Vector* y = 0);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_GRG_H
