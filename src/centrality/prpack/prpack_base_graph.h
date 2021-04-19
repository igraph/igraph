#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csc.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <cstdio>
#include <utility>

namespace prpack {

    class prpack_base_graph {
        private:
            // helper methods
            void initialize();
            void read_smat(std::FILE* f, const bool weighted);
            void read_edges(std::FILE* f);
            void read_ascii(std::FILE* f);
        public:
            // instance variables
            int num_vs;
            int num_es;
            int num_self_es;
            int* heads;
            int* tails;
            double* vals;
            // constructors
            prpack_base_graph();    // only to support inheritance
            prpack_base_graph(const prpack_csc* g);
            prpack_base_graph(const prpack_int64_csc* g);
            prpack_base_graph(const prpack_csr* g);
            prpack_base_graph(const prpack_edge_list* g);
            prpack_base_graph(const char* filename, const char* format, const bool weighted);
            prpack_base_graph(int nverts, int nedges, std::pair<int,int>* edges);
            // destructor
            ~prpack_base_graph();
            // operations
            void normalize_weights();
    };

}

#endif
