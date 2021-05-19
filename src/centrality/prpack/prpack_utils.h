#ifndef PRPACK_UTILS
#define PRPACK_UTILS
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#include <string>

// Computes the time taken to do X and stores it in T.
#define TIME(T, X)                  \
    (T) = prpack_utils::get_time(); \
    (X);                            \
    (T) = prpack_utils::get_time() - (T)

// Computes S += A using C as a carry-over.
// This is a macro over a function as it is faster this way.
#define COMPENSATED_SUM(S, A, C)                        \
    double compensated_sum_y = (A) - (C);               \
    double compensated_sum_t = (S) + compensated_sum_y; \
    (C) = compensated_sum_t - (S) - compensated_sum_y;  \
    (S) = compensated_sum_t

namespace prpack {

    class prpack_utils {
        public:
            static double get_time();
            static void validate(const bool condition, const std::string& msg);
            static double* permute(const int length, const double* a, const int* coding);
    };

}

#endif
