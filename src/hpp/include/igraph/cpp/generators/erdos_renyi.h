/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_ERDOS_RENYI_H
#define IGRAPHPP_GENERATORS_ERDOS_RENYI_H

#include <memory>

namespace igraph {

class Graph;

/// Generates an Erdos-Renyi random network with the G(n, m) model
std::auto_ptr<Graph> erdos_renyi_game_gnm(integer_t n, integer_t m,
        bool directed = false, bool loops = false);

/// Generates an Erdos-Renyi random network with the G(n, p) model
std::auto_ptr<Graph> erdos_renyi_game_gnp(integer_t n, real_t p,
        bool directed = false, bool loops = false);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_ERDOS_RENYI_H
