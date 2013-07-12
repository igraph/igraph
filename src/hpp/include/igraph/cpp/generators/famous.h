/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_FAMOUS_H
#define IGRAPHPP_GENERATORS_FAMOUS_H

#include <memory>

namespace igraph {

class Graph;

/// Generates a famous graph given its name
std::auto_ptr<Graph> famous(const std::string& name);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_FAMOUS_H
