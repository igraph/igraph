/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/graph.h>
#include <igraph/cpp/generators/famous.h>

namespace igraph {

std::auto_ptr<Graph> famous(const std::string& name) {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_famous(result.get(), name.c_str()));
    return std::auto_ptr<Graph>(new Graph(result));
}

}         // end of namespaces

