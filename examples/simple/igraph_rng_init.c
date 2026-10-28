#include <igraph.h>
#include <stdio.h>

int main(void) {

    /* Initialize the library. */
    igraph_setup();

    /* 1. Create a new random number generator. */
    igraph_rng_t rng;

    /* Initialize the RNG using the MT19937 algorithm. */
    igraph_rng_init(&rng, &igraph_rngtype_mt19937);

    /* 2. Generate some random numbers with it. */
    igraph_int_t roll;

    for (int i = 0; i < 5; i++) {
        roll = igraph_rng_get_integer(&rng, 1, 6);
        printf("Die roll: %" IGRAPH_PRId "\n", roll);
    }

    /* 3. Make it the default random number generator. */
    /* Save the previous default RNG so we can restore it later. */
    igraph_rng_t *old_rng = igraph_rng_set_default(&rng);

    /* Get random numbers using the current default RNG. */
    for (int i = 0; i < 5; i++) {
        roll = igraph_rng_get_integer(igraph_rng_default(), 1, 6);
        printf("Die roll (default RNG): %" IGRAPH_PRId "\n", roll);
    }

    /* Clean up. */
    /* Restore the previous default RNG before destroying our RNG. */
    igraph_rng_set_default(old_rng);
    igraph_rng_destroy(&rng);

    return 0;
}
