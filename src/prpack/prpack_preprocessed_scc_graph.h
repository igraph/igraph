#ifndef PRPACK_PREPROCESSED_SCC_GRAPH
#define PRPACK_PREPROCESSED_SCC_GRAPH
#include "prpack_preprocessed_graph.h"
#include "prpack_base_graph.h"

namespace prpack {

    // Pre-processed graph class
    class prpack_preprocessed_scc_graph : public prpack_preprocessed_graph {
        private:
            // helper methods
            void initialize();
            void initialize_weighted(const prpack_base_graph* bg);
            void initialize_unweighted(const prpack_base_graph* bg);
        public:
            // instance variables
            int num_es_inside;
            int* heads_inside;
            int* tails_inside;
            double* vals_inside;
            int num_es_outside;
            int* heads_outside;
            int* tails_outside;
            double* vals_outside;
            double* ii;
            double* num_outlinks;
            int num_comps;
            int* divisions;
            int* encoding;
            int* decoding;
            // constructors
            prpack_preprocessed_scc_graph(const prpack_base_graph* bg);
            // destructor
            ~prpack_preprocessed_scc_graph();
    };

}

#endif
