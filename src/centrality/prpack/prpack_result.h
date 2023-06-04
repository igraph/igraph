#ifndef PRPACK_RESULT
#define PRPACK_RESULT

#include <string>
#include <cstdint>

namespace prpack {

    // Result class.
    class prpack_result {
        public:
            // instance variables
            int num_vs;
            int num_es;
            double* x;
            double read_time;
            double preprocess_time;
            double compute_time;
            int64_t num_es_touched;
            std::string method;
            int converged;
            // constructor
            prpack_result();
            // destructor
            ~prpack_result();
    };

}

#endif
