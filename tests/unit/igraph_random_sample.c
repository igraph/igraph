/* igraph library.
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

void check(igraph_int_t l, igraph_int_t h, igraph_int_t s) {
    igraph_vector_int_t v;

    igraph_vector_int_init(&v, 0);
    igraph_random_sample(&v, l, h, s);

    IGRAPH_ASSERT(igraph_vector_int_size(&v) == s);
    IGRAPH_ASSERT(VECTOR(v)[0] >= l);
    IGRAPH_ASSERT(VECTOR(v)[igraph_vector_int_size(&v) - 1] <= h);

    for (igraph_int_t i = 0; i < s - 1; i++) {
        IGRAPH_ASSERT(VECTOR(v)[i] < VECTOR(v)[i + 1]);
    }

    igraph_vector_int_destroy(&v);
}

int main(void) {
    igraph_vector_int_t v;

    igraph_rng_seed(igraph_rng_default(), 42);

    check(-100, 100, 10);
    check(-10000, 10000, 100);
    check(0, 0, 1);
    check(100, 110, 11);

    igraph_vector_int_init(&v, 0);

    igraph_random_sample(&v, 100, 1000, 0);
    IGRAPH_ASSERT(igraph_vector_int_size(&v) == 0);

    /* test parameters */
    /*----------low----high----length*/
    /* lower limit is greater than upper limit */
    CHECK_ERROR(igraph_random_sample(&v, 300, 200, 10), IGRAPH_EINVAL);
    /* sample size is greater than size of candidate pool */
    CHECK_ERROR(igraph_random_sample(&v, 200, 300, 500), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
