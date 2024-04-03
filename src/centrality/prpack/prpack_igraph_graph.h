#ifndef PRPACK_IGRAPH_GRAPH
#define PRPACK_IGRAPH_GRAPH

#ifdef PRPACK_IGRAPH_SUPPORT

#include "prpack_base_graph.h"

#include "igraph_datatype.h"
#include "igraph_vector.h"

namespace prpack {

    class prpack_igraph_graph : public prpack_base_graph {
    public:
        // constructors
        prpack_igraph_graph() { }

        // We use a separate function to carry out the actual construction of the graph.
        // The base class constructor sets the heads/tails/vals arrays to NULL,
        // so these can safely be delete'ed by the destructor when
        // convert_from_igraph() fails.
        igraph_error_t convert_from_igraph(const igraph_t *g,
                                           const igraph_vector_t *weights,
                                           bool directed = true);
    };

}

// PRPACK_IGRAPH_SUPPORT
#endif

// PRPACK_IGRAPH_GRAPH
#endif
