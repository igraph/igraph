
#include <igraph.h>

#include "bench.h"

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

void bench(int m, int n, int rep) {
    igraph_matrix_t mat;

    igraph_matrix_init(&mat, m, n);
    RNG_BEGIN();
    for (igraph_integer_t j=0; j < n; j++) {
        for (igraph_integer_t i=0; i < m; i++) {
            MATRIX(mat, i, j) = RNG_UNIF(-1, 1);
        }
    }
    RNG_END();

    char name[200];
    sprintf(name, "Transpose %5d x %5d, %dx", m, n, rep);
    BENCH( name, REPEAT(igraph_matrix_transpose(&mat), rep) );

    igraph_matrix_destroy(&mat);
}

int main(void) {
    bench(30, 30, 100000);
    bench(100, 100, 10000);
    bench(1000, 1000, 100);
    bench(1024, 1024, 100); /* naive implementation has bad cache behaviour with power of 2 sizes */
    bench(1023, 1025, 100); /* non-symmetric */
    bench(3000, 3000, 10);
    /* skinny non-symmetric: */
    bench(100, 10000, 100);
    bench(10000, 100, 100);

    return 0;
}
