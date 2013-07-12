/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/cpp/graph.h>
#include <igraph/cpp/generators/full.h>

namespace igraph {

std::auto_ptr<Graph> full(integer_t nodes, bool directed, bool loops) {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_full(result.get(), nodes, directed, loops));
    return std::auto_ptr<Graph>(new Graph(result));
}

}         // end of namespaces

