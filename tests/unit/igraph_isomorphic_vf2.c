/* IGraph library.
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

int main(void) {
    igraph_t ring1, ring2;
    igraph_vector_int_t color1, color2;
    igraph_vector_int_t perm;
    igraph_bool_t iso;
    igraph_integer_t count;
    igraph_integer_t i;

    igraph_rng_seed(igraph_rng_default(), 12345);

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/1);
    igraph_vector_int_init_range(&perm, 0, igraph_vcount(&ring1));
    igraph_vector_int_shuffle(&perm);
    igraph_permute_vertices(&ring1, &ring2, &perm);

    /* Without colors */
    igraph_isomorphic(&ring1, &ring2, &iso);
    if (!iso) {
        fprintf(stderr, "Without color failed.\n");
        return 1;
    }

    /* Without colors, number of isomorphisms */
    igraph_count_isomorphisms_vf2(&ring1, &ring2, 0, 0, 0, 0, &count, 0, 0, 0);
    if (count != 200) {
        fprintf(stderr, "Count without colors failed, expected 200, got %" IGRAPH_PRId ".\n", count);
        return 2;
    }

    /* Separate colors for each vertex */
    igraph_vector_int_init(&color1, igraph_vcount(&ring1));
    igraph_vector_int_init(&color2, igraph_vcount(&ring2));
    for (i = 0; i < igraph_vector_int_size(&color1); i++) {
        VECTOR(color1)[i] = VECTOR(color2)[VECTOR(perm)[i]] = i;
    }
    igraph_count_isomorphisms_vf2(&ring1, &ring2, &color1, &color2, 0, 0, &count, 0, 0, 0);
    if (count != 1) {
        fprintf(stderr, "Count with separate colors failed, expected 1, got %" IGRAPH_PRId ".\n", count);
        return 5;
    }

    /* Try a negative result */
    igraph_vector_int_fill(&color1, 0);
    igraph_vector_int_fill(&color2, 0);
    VECTOR(color1)[0] = 1;
    igraph_isomorphic_vf2(&ring1, &ring2, &color1, &color2, 0, 0, &iso, 0, 0, 0, 0, 0);
    if (iso) {
        fprintf(stderr, "Negative test failed.\n");
        return 6;
    }

    /* Another negative, same color distribution, different topology */
    igraph_vector_int_fill(&color1, 0);
    igraph_vector_int_fill(&color2, 0);
    VECTOR(color1)[0] = 1;
    VECTOR(color1)[1] = 1;
    VECTOR(color2)[0] = 1;
    VECTOR(color2)[(VECTOR(perm)[1] + 1) % igraph_vcount(&ring2)] = 1;
    igraph_isomorphic_vf2(&ring1, &ring2, &color1, &color2, 0, 0, &iso, 0, 0, 0, 0, 0);
    if (iso) {
        fprintf(stderr, "Second negative test failed.\n");
        return 7;
    }

    igraph_vector_int_destroy(&color1);
    igraph_vector_int_destroy(&color2);

    igraph_vector_int_destroy(&perm);
    igraph_destroy(&ring2);
    igraph_destroy(&ring1);

    /* ---------------------------------------------------------------- */
    /* SUBGRAPH ISOMORPHISM                                             */
    /* ---------------------------------------------------------------- */

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/0);
    igraph_ring(&ring2, 80, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/0);

    /* One color */
    igraph_vector_int_init(&color1, igraph_vcount(&ring1));
    igraph_vector_int_init(&color2, igraph_vcount(&ring2));
    igraph_count_subisomorphisms_vf2(&ring1, &ring2, &color1, &color2, 0, 0,
                                     &count, 0, 0, 0);
    if (count != 42) {
        fprintf(stderr, "Count with one color failed, expected 42, got %" IGRAPH_PRId ".\n", count);
        return 31;
    }

    igraph_vector_int_destroy(&color1);
    igraph_vector_int_destroy(&color2);

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);

    /* ---------------------------------------------------------------- */
    /* EDGE COLORING, GRAPH ISOMORPHISM                                 */
    /* ---------------------------------------------------------------- */

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_vector_int_init_range(&perm, 0, igraph_ecount(&ring1));
    igraph_vector_int_shuffle(&perm);
    igraph_permute_vertices(&ring1, &ring2, &perm);
    igraph_vector_int_destroy(&perm);

    /* Everything has the same color */
    igraph_vector_int_init(&color1, igraph_ecount(&ring1));
    igraph_vector_int_init(&color2, igraph_ecount(&ring2));
    igraph_isomorphic_vf2(&ring1, &ring2, 0, 0, &color1, &color2, &iso, 0, 0, 0, 0, 0);
    if (!iso) {
        fprintf(stderr, "Single edge-color failed.\n");
        return 41;
    }

    /* Two colors, just counting */
    for (i = 0; i < igraph_vector_int_size(&color1); i += 2) {
        VECTOR(color1)[i]   = VECTOR(color2)[i] = 0;
        VECTOR(color1)[i + 1] = VECTOR(color2)[i] = 1;
    }
    igraph_count_isomorphisms_vf2(&ring1, &ring2, 0, 0, &color1, &color2, &count, 0, 0, 0);
    if (count != 100) {
        fprintf(stderr, "Count with two edge colors failed, expected 100, got %" IGRAPH_PRId ".\n", count);
        return 42;
    }

    /* Separate colors for each edge */
    for (i = 0; i < igraph_vector_int_size(&color1); i++) {
        VECTOR(color1)[i]   = VECTOR(color2)[i] = i;
    }
    igraph_count_isomorphisms_vf2(&ring1, &ring2, 0, 0, &color1, &color2, &count, 0, 0, 0);
    if (count != 1) {
        fprintf(stderr, "Count with separate edge colors failed, expected 1, got %" IGRAPH_PRId ".\n", count);
        return 43;
    }

    /* Try a negative result */
    igraph_vector_int_fill(&color1, 0);
    igraph_vector_int_fill(&color2, 0);
    VECTOR(color1)[0] = 1;
    igraph_isomorphic_vf2(&ring1, &ring2, 0, 0, &color1, &color2, &iso, 0, 0, 0, 0, 0);
    if (iso) {
        fprintf(stderr, "Negative edge test failed.\n");
        return 44;
    }

    /* Another negative, same color distribution, different topology */
    igraph_vector_int_fill(&color1, 0);
    igraph_vector_int_fill(&color2, 0);
    VECTOR(color1)[0] = 1;
    VECTOR(color1)[1] = 1;
    VECTOR(color2)[0] = 1;
    VECTOR(color2)[2] = 1;
    igraph_isomorphic_vf2(&ring1, &ring2, 0, 0, &color1, &color2, &iso, 0, 0, 0, 0, 0);
    if (iso) {
        fprintf(stderr, "Second negative edge test failed.\n");
        return 45;
    }

    igraph_vector_int_destroy(&color1);
    igraph_vector_int_destroy(&color2);

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);

    /* ---------------------------------------------------------------- */
    /* EDGE COLORED SUBGRAPH ISOMORPHISM                                */
    /* ---------------------------------------------------------------- */

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/0);
    igraph_ring(&ring2, 80, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/0);

    /* One color */
    igraph_vector_int_init(&color1, igraph_ecount(&ring1));
    igraph_vector_int_init(&color2, igraph_ecount(&ring2));
    igraph_count_subisomorphisms_vf2(&ring1, &ring2, 0, 0, &color1, &color2,
                                     &count, 0, 0, 0);
    if (count != 42) {
        fprintf(stderr, "Count with one edge color failed, expected 42, got %" IGRAPH_PRId ".\n", count);
        return 51;
    }

    /* Two colors */
    for (i = 0; i < igraph_vector_int_size(&color1) - 1; i += 2) {
        VECTOR(color1)[i]   = 0;
        VECTOR(color1)[i + 1] = 1;
    }
    for (i = 0; i < igraph_vector_int_size(&color2) - 1; i += 2) {
        VECTOR(color2)[i]   = 0;
        VECTOR(color2)[i + 1] = 1;
    }
    igraph_count_subisomorphisms_vf2(&ring1, &ring2, 0, 0, &color1, &color2,
                                     &count, 0, 0, 0);
    if (count != 22) {
        fprintf(stderr, "Count with two edge colors failed, expected 22, got %" IGRAPH_PRId ".\n", count);
        return 52;
    }

    igraph_vector_int_destroy(&color1);
    igraph_vector_int_destroy(&color2);

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);

    VERIFY_FINALLY_STACK();

    return 0;
}
