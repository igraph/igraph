/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/igraph_components.h>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/analysis/components.h>

namespace igraph {

void clusters(const Graph& graph, Vector* membership,
        Vector* csize, integer_t* no, Connectedness connectedness) {
    igraph_vector_t* c_membership;
    igraph_vector_t* c_csize;

    c_membership = (membership == 0) ? 0 : membership->c_vector();
    c_csize = (csize == 0) ? 0 : csize->c_vector();

    IGRAPH_TRY(igraph_clusters(graph.c_graph(), c_membership, c_csize,
                no, connectedness));
}

}         // end of namespaces
