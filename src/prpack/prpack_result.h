#ifndef PRPACK_RESULT
#define PRPACK_RESULT

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
            long num_es_touched;
            const char* method;
            int converged;
            // constructor
            prpack_result();
            // destructor
            ~prpack_result();
    };

};

#endif
