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
        explicit prpack_igraph_graph(const igraph_t *g,
                                     const igraph_vector_t *weights = 0,
                                     bool directed = true);
    };

}

// PRPACK_IGRAPH_SUPPORT
#endif

// PRPACK_IGRAPH_GRAPH
#endif
