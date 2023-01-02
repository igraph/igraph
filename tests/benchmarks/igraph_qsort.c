#include <igraph.h>

#include "bench.h"

/* This program benchmarks igraph_qsort() indirectly through vector_sort() */

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

int main(void) {
    igraph_vector_int_t vec;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();


    igraph_vector_int_init(&vec, 0);

    RNG_BEGIN();

#define N 10000000

    igraph_vector_int_resize(&vec, N);
    for (igraph_integer_t i=0; i < N; i++) {
        VECTOR(vec)[i] = RNG_INTEGER(0, N-1);
    }
    BENCH("Sort vector of length " TOSTRING(N), igraph_vector_int_sort(&vec));

#undef N
#define N 1000000

    igraph_vector_int_resize(&vec, N);
    for (igraph_integer_t i=0; i < N; i++) {
        VECTOR(vec)[i] = RNG_INTEGER(0, N-1);
    }
    BENCH("Sort vector of length " TOSTRING(N), igraph_vector_int_sort(&vec));

#undef N
#define N 100000

    igraph_vector_int_resize(&vec, N);
    for (igraph_integer_t i=0; i < N; i++) {
        VECTOR(vec)[i] = RNG_INTEGER(0, N-1);
    }
    BENCH("Sort vector of length " TOSTRING(N), igraph_vector_int_sort(&vec));

    RNG_END();

    igraph_vector_int_destroy(&vec);

    return 0;
}
