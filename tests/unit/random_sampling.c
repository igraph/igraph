/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

/* These is merely a smoke test for various random samplers.
 * It does not verify the correctness of the result, except
 * for some special edge cases. */

void sample() {
    igraph_integer_t i;
    igraph_real_t x;

    i = RNG_INTEGER(-100, 100);
    printf("integer: %" IGRAPH_PRId "\n", i);
    IGRAPH_ASSERT(-100 <= i && i <= 100);

    i = RNG_INTEGER(2, 2);
    printf("integer: %" IGRAPH_PRId "\n", i);
    IGRAPH_ASSERT(i == 2);

    i = RNG_INTEGER(-5, -5);
    printf("integer: %" IGRAPH_PRId "\n", i);
    IGRAPH_ASSERT(i == -5);

    x = RNG_UNIF(-100, 100);
    printf("unif: %g\n", x);
    IGRAPH_ASSERT(-100 <= x && x < 100);

    x = RNG_UNIF(3, 3);
    printf("unif: %g\n", x);
    IGRAPH_ASSERT(x == 3);

    x = RNG_UNIF(-4.5, -4.5);
    printf("unif: %g\n", x);
    IGRAPH_ASSERT(x == -4.5);

    x = RNG_UNIF01();
    printf("unif01: %g\n", x);
    IGRAPH_ASSERT(0 <= x);
    IGRAPH_ASSERT(x < 1);

    x = RNG_BINOM(10, 0.3);
    printf("binom: %g\n", x);
    IGRAPH_ASSERT(0 <= x && x <= 10);

    x = RNG_BINOM(5, 0);
    IGRAPH_ASSERT(x == 0);

    x = RNG_BINOM(5, 1);
    IGRAPH_ASSERT(x == 5);

    x = RNG_GEOM(0.2);
    printf("geom: %g\n", x);
    IGRAPH_ASSERT(0 <= x);

    x = RNG_GEOM(1);
    IGRAPH_ASSERT(x == 0);

    x = RNG_GEOM(0);
    IGRAPH_ASSERT(igraph_is_nan(x));

    x = RNG_NORMAL(0, 1);
    printf("normal: %g\n", x);
    IGRAPH_ASSERT(IGRAPH_FINITE(x));

    x = RNG_EXP(1);
    printf("exp: %g\n", x);
    IGRAPH_ASSERT(0 <= x);

    x = RNG_EXP(0);
    IGRAPH_ASSERT(igraph_is_nan(x));

    x = RNG_GAMMA(1, 1);
    printf("gamma: %g\n", x);
    IGRAPH_ASSERT(0 <= x);

    /* Note: Some systems would reject 0 as a parameter value instead of returning 0 */
    x = RNG_GAMMA(0, 1);
    IGRAPH_ASSERT(x == 0);

    x = RNG_GAMMA(1, 0);
    IGRAPH_ASSERT(x == 0);

    x = RNG_POIS(10);
    printf("poisson: %g\n", x);
    IGRAPH_ASSERT(0 <= x && IGRAPH_FINITE(x));

    x = RNG_POIS(0);
    IGRAPH_ASSERT(x == 0);

    x = RNG_POIS(-1);
    IGRAPH_ASSERT(igraph_is_nan(x));
}

void test_and_destroy(igraph_rng_type_t *rng_type) {
    igraph_rng_t *def = igraph_rng_default();
    igraph_rng_t rng;

    igraph_error_handler_t *oldhandler = igraph_set_error_handler(&igraph_error_handler_printignore);
    igraph_error_t err = igraph_rng_init(&rng, rng_type);
    switch (err) {
    case IGRAPH_SUCCESS:
        break;
    case IGRAPH_UNIMPLEMENTED:
        return;
    default:
        IGRAPH_FATAL("Error while initializing RNG.");
    }
    igraph_set_error_handler(oldhandler);

    printf("\n%s\n\n", igraph_rng_name(&rng));

    igraph_rng_set_default(&rng);
    igraph_rng_seed(igraph_rng_default(), 137);

    sample();

    igraph_rng_set_default(def);
    igraph_rng_destroy(&rng);
}

int main() {
    igraph_rng_type_t rng_types[] = {
        igraph_rngtype_glibc2,
        igraph_rngtype_mt19937,
        igraph_rngtype_pcg32,
        igraph_rngtype_pcg64
    };
    igraph_integer_t i;

    printf("Default\n\n");
    igraph_rng_seed(igraph_rng_default(), 709);

    sample();

    for (i = 0; i < sizeof(rng_types) / sizeof(rng_types[0]); i++) {
        test_and_destroy(&rng_types[i]);
    }

    return 0;
}
