
#include <igraph.h>

/*
 * This test serves to ensure that the same sequence of random numbers are generated for the
 * same seed on all platforms (different operating systems and 32- or 64-bit systems).
 */

int main() {
    int i;
    igraph_rng_seed(igraph_rng_default(), 137);

    for (i = 0; i < 32; ++i) {
        printf("%ld\n", RNG_INTEGER(0, 100));
    }

    for (i = 0; i < 32; ++i) {
        printf("%g\n", RNG_UNIF(0, 1e-6));
    }
    return 0;
}
