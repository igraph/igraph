#ifndef PRPACK_EDGE_LIST
#define PRPACK_EDGE_LIST

namespace prpack {

    class prpack_edge_list {
        public:
            int num_vs;
            int num_es;
            int* heads;
            int* tails;
    };

}

#endif
