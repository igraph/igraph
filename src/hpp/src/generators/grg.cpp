/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/igraph_games.h>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/generators/grg.h>

namespace igraph {

std::auto_ptr<Graph> grg_game(integer_t nodes, real_t radius, bool torus,
        Vector* x, Vector* y) {
    std::auto_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_grg_game(result.get(), nodes, radius, torus,
                x ? x->c_vector() : 0, y ? y->c_vector() : 0));
    return std::auto_ptr<Graph>(new Graph(result));
}

}         // end of namespaces

