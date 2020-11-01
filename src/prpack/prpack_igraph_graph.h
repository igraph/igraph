#ifndef PRPACK_IGRAPH_GRAPH
#define PRPACK_IGRAPH_GRAPH

#ifdef PRPACK_IGRAPH_SUPPORT

#include "igraph_interface.h"
#include "prpack_base_graph.h"

namespace prpack {

    class prpack_igraph_graph : public prpack_base_graph {

        public:
            // constructors
            explicit prpack_igraph_graph(const igraph_t* g,
					const igraph_vector_t* weights = 0,
					igraph_bool_t directed = true);
    };

}

// PRPACK_IGRAPH_SUPPORT 
#endif 

// PRPACK_IGRAPH_GRAPH
#endif
