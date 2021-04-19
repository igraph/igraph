#ifndef PRPACK_IGRAPH_GRAPH
#define PRPACK_IGRAPH_GRAPH

#ifdef PRPACK_IGRAPH_SUPPORT

#include "prpack_base_graph.h"

struct igraph_s;
struct igraph_vector_t;

namespace prpack {

    class prpack_igraph_graph : public prpack_base_graph {

        public:
            // constructors
            explicit prpack_igraph_graph(const struct igraph_s* g,
					const struct igraph_vector_t* weights = 0,
					bool directed = true);
    };

}

// PRPACK_IGRAPH_SUPPORT
#endif

// PRPACK_IGRAPH_GRAPH
#endif
