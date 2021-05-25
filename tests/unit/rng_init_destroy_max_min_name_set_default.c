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
#include "test_utilities.inc"

void test_and_destroy(igraph_rng_type_t *rng_type, igraph_rng_t *rng_def) {
    int i;
    igraph_rng_t rng;

    IGRAPH_ASSERT(igraph_rng_init(&rng, rng_type) == IGRAPH_SUCCESS);
    printf("rng name: %s\n", igraph_rng_name(&rng));

    igraph_rng_seed(&rng, 42);
    for (i = 0; i < 5; i++) {
        printf("%ld\n", igraph_rng_get_integer(&rng, 0, 100));
    }
    printf("\n");

    igraph_rng_set_default(&rng);
    igraph_rng_seed(igraph_rng_default(), 42);
    for (i = 0; i < 5; i++) {
        printf("%ld\n", igraph_rng_get_integer(igraph_rng_default(), 0, 100));
    }
    printf("\n");

    IGRAPH_ASSERT(igraph_rng_max(&rng) >= 32767);
    igraph_rng_set_default(rng_def);
    igraph_rng_destroy(&rng);
}

int main() {
    int i;
    igraph_rng_type_t rng_types[3] = {igraph_rngtype_glibc2, igraph_rngtype_mt19937, igraph_rngtype_rand};
    igraph_rng_t rng_def;

    IGRAPH_ASSERT(igraph_rng_init(&rng_def, &igraph_rngtype_glibc2) == IGRAPH_SUCCESS);

    for (i = 0; i < 3; i++) {
        test_and_destroy(&rng_types[i], &rng_def);
    }

    igraph_rng_destroy(&rng_def);

    VERIFY_FINALLY_STACK();
    return 0;
}
