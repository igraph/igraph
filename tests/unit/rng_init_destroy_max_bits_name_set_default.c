/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void test_and_destroy(igraph_rng_type_t *rng_type, igraph_rng_t *rng_def) {
    int i;
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

    printf("rng name: %s\n", igraph_rng_name(&rng));

    igraph_rng_seed(&rng, 42);
    for (i = 0; i < 5; i++) {
        printf("%" IGRAPH_PRId "\n", igraph_rng_get_integer(&rng, 0, 100));
    }
    printf("\n");

    igraph_rng_set_default(&rng);
    igraph_rng_seed(igraph_rng_default(), 42);
    for (i = 0; i < 5; i++) {
        printf("%" IGRAPH_PRId "\n", igraph_rng_get_integer(igraph_rng_default(), 0, 100));
    }
    printf("\n");

    IGRAPH_ASSERT(igraph_rng_max(&rng) >= 0x7fffffff);
    IGRAPH_ASSERT(igraph_rng_bits(&rng) >= 31);

    igraph_rng_set_default(rng_def);
    igraph_rng_destroy(&rng);
}

void test_and_destroy_with_expected_values(
    igraph_rng_type_t *rng_type, const igraph_vector_int_t *expected
) {
    int i;
    igraph_error_t retval;
    igraph_rng_t rng;

    igraph_set_error_handler(igraph_error_handler_ignore);
    retval = igraph_rng_init(&rng, rng_type);
    igraph_set_error_handler(igraph_error_handler_abort);

    if (retval == IGRAPH_UNIMPLEMENTED) {
        /* not supported in the current setup, this is OK */
        return;
    }

    igraph_rng_seed(&rng, 42);
    for (i = 0; i < 5; i++) {
        IGRAPH_ASSERT(VECTOR(*expected)[i] == igraph_rng_get_integer(&rng, 0, 100));
    }

    igraph_rng_seed(&rng, 42);
    for (i = 0; i < 5; i++) {
        IGRAPH_ASSERT(VECTOR(*expected)[i] == igraph_rng_get_integer(&rng, 0, 100));
    }

    IGRAPH_ASSERT(igraph_rng_max(&rng) >= 0x7fffffff);
    IGRAPH_ASSERT(igraph_rng_bits(&rng) >= 31);

    igraph_rng_destroy(&rng);
}

void test_mandatory_rngtypes(void) {
    int i;
    igraph_rng_type_t rng_types[] = {
        igraph_rngtype_glibc2,
        igraph_rngtype_mt19937,
        igraph_rngtype_pcg32,
    };
    const int NUM_RNG_TYPES = sizeof(rng_types) / sizeof(rng_types[0]);
    igraph_rng_t rng_def;

    IGRAPH_ASSERT(igraph_rng_init(&rng_def, &igraph_rngtype_glibc2) == IGRAPH_SUCCESS);

    for (i = 0; i < NUM_RNG_TYPES; i++) {
        test_and_destroy(&rng_types[i], &rng_def);
    }

    igraph_rng_destroy(&rng_def);

    VERIFY_FINALLY_STACK();
}

void test_optional_rngtypes(void) {
    igraph_rng_type_t rng_types[] = {
        igraph_rngtype_pcg64
    };
    igraph_vector_int_t expected;
    igraph_integer_t expected_values[1][5] = {
        { 61, 38, 73, 84, 67 }
    };
    int i;

    for (i = 0; i < 1; i++) {
        igraph_vector_int_view(&expected, expected_values[i], 5);
        test_and_destroy_with_expected_values(&rng_types[i], &expected);
    }
}

int main(void) {
    test_mandatory_rngtypes();
    test_optional_rngtypes();
    return 0;
}
