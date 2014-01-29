#ifndef PRPACK_CSC
#define PRPACK_CSC

#include <stdint.h>

namespace prpack {

    class prpack_csc {
        public:
            int num_vs;
            int num_es;
            int* heads;
            int* tails;
    };

    class prpack_int64_csc {
        public:
            int64_t num_vs;
            int64_t num_es;
            int64_t* heads;
            int64_t* tails;
    };
};

#endif
