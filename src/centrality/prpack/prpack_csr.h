#ifndef PRPACK_CSR
#define PRPACK_CSR

namespace prpack {

    class prpack_csr {
        public:
            int num_vs;
            int num_es;
            int* heads;
            int* tails;
    };

}

#endif
