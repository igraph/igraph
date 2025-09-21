/*
   igraph library.
   Copyright (C) 2015-2022  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>
#include <string.h>

#include "test_utilities.h"

void random_permutation(igraph_vector_int_t *vec) {
    /* We just do size(vec) * 2 swaps */
    igraph_int_t one, two, tmp, i, n = igraph_vector_int_size(vec);
    for (i = 0; i < 2 * n; i++) {
        one = RNG_INTEGER(0, n - 1);
        two = RNG_INTEGER(0, n - 1);
        tmp = VECTOR(*vec)[one];
        VECTOR(*vec)[one] = VECTOR(*vec)[two];
        VECTOR(*vec)[two] = tmp;
    }
}


void test3(void) {
    igraph_int_t i, j;
    igraph_graph_list_t graphs3;
    igraph_t g;

    // Verify that no two 3-vertex graphs of distinct isoclasses are considered isomorphic by Bliss or VF2.

    igraph_graph_list_init(&graphs3, 0);

    for (i = 0; i < 16; i++) {
        igraph_isoclass_create(&g, 3, i, /* directed = */ 1);
        igraph_graph_list_push_back(&graphs3, &g);
    }

    for (i = 0; i < 16; i++)
        for (j = i + 1; j < 16; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_bliss(
                igraph_graph_list_get_ptr(&graphs3, i),
                igraph_graph_list_get_ptr(&graphs3, j),
                NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL
            );
            if (iso) {
                printf("Bliss failure, 3 vertex directed graphs of isoclass %" IGRAPH_PRId " and %" IGRAPH_PRId " are not isomorphic. Bliss reports otherwise.\n", i, j);
            }
        }

    for (i = 0; i < 16; i++)
        for (j = i + 1; j < 16; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_vf2(
                igraph_graph_list_get_ptr(&graphs3, i),
                igraph_graph_list_get_ptr(&graphs3, j),
                NULL, NULL, NULL, NULL, &iso, NULL, NULL, NULL, NULL, NULL
            );
            if (iso) {
                printf("VF2 failure, 3 vertex directed graphs of isoclass %" IGRAPH_PRId " and %" IGRAPH_PRId " are not isomorphic. VF2 reports otherwise.\n", i, j);
            }
        }

    igraph_graph_list_destroy(&graphs3);
}


void test4(void) {
    igraph_int_t i, j;
    igraph_graph_list_t graphs4;
    igraph_t g;

    // Verify that no two 4-vertex graphs of distinct isoclasses are considered isomorphic by Bliss or VF2.

    igraph_graph_list_init(&graphs4, 0);

    for (i = 0; i < 218; i++) {
        igraph_isoclass_create(&g, 4, i, /* directed = */ 1);
        igraph_graph_list_push_back(&graphs4, &g);
    }

    for (i = 0; i < 218; i++)
        for (j = i + 1; j < 218; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_bliss(
                igraph_graph_list_get_ptr(&graphs4, i),
                igraph_graph_list_get_ptr(&graphs4, j),
                NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL
            );
            if (iso) {
                printf("Bliss failure, 4 vertex directed graphs of isoclass %" IGRAPH_PRId " and %" IGRAPH_PRId " are not isomorphic. Bliss reports otherwise.\n", i, j);
            }
        }

    for (i = 0; i < 218; i++)
        for (j = i + 1; j < 218; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_vf2(
                igraph_graph_list_get_ptr(&graphs4, i),
                igraph_graph_list_get_ptr(&graphs4, j),
                NULL, NULL, NULL, NULL, &iso, NULL, NULL, NULL, NULL, NULL
            );
            if (iso) {
                printf("VF2 failure, 4 vertex directed graphs of isoclass %" IGRAPH_PRId " and %" IGRAPH_PRId " are not isomorphic. VF2 reports otherwise.\n", i, j);
            }
        }

    igraph_graph_list_destroy(&graphs4);
}


void test_bliss(void) {
    igraph_t ring1, ring2, directed_ring;
    igraph_vector_int_t perm;
    igraph_bool_t iso;
    igraph_bliss_info_t info;
    igraph_vector_int_t color;
    igraph_vector_int_list_t generators;
    igraph_real_t num_automorphisms;

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/1);
    igraph_vector_int_init_range(&perm, 0, igraph_vcount(&ring1));
    random_permutation(&perm);
    igraph_permute_vertices(&ring1, &ring2, &perm);

    igraph_ring(&directed_ring, 100, /* directed= */ 1, /* mutual = */0, /* circular = */1);

    igraph_vector_int_list_init(&generators, 0);

    igraph_isomorphic_bliss(&ring1, &ring2, NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
    if (! iso) {
        printf("Bliss failed on ring isomorphism.\n");
    }

    igraph_count_automorphisms(&ring1, NULL, &num_automorphisms);
    if (num_automorphisms != 200) {
        printf("Biss automorphism count failed: ring1.\n");
    }

    igraph_count_automorphisms(&ring2, NULL, &num_automorphisms);
    if (num_automorphisms != 200) {
        printf("Biss automorphism count failed: ring2.\n");
    }

    igraph_count_automorphisms(&directed_ring, NULL, &num_automorphisms);
    if (num_automorphisms != 100) {
        printf("Biss automorphism count failed: directed_ring.\n");
    }

    // The following test is included so there is at least one call to igraph_automorphism_group_bliss
    // in the test suite. However, the generator set returned may depend on the splitting
    // heursitics as well as on the Bliss version. If the test fails, please verify manually
    // that the generating set is valid. For a undirected cycle graph like ring2, there should
    // be two generators: a cyclic permutation and a reversal of the vertex order.
    igraph_automorphism_group_bliss(&ring2, NULL, &generators, IGRAPH_BLISS_F, NULL);
    if (igraph_vector_int_list_size(&generators) != 2)
        printf("Bliss automorphism generators may have failed with ring2. "
               "Please verify the generators manually. "
               "Note that the generator set is not guaranteed to be minimal.\n");
    igraph_vector_int_list_clear(&generators);

    // For a directed ring, the only generator should be a cyclic permutation.
    igraph_automorphism_group_bliss(&directed_ring, NULL, &generators, IGRAPH_BLISS_F, NULL);
    if (igraph_vector_int_list_size(&generators) != 1)
        printf("Bliss automorphism generators may have failed with directed_ring. "
               "Please verify the generators manually. "
               "Note that the generator set is not guaranteed to be minimal.\n");
    igraph_vector_int_list_clear(&generators);

    igraph_vector_int_init_range(&color, 0, igraph_vcount(&ring1));

    igraph_count_automorphisms_bliss(&ring1, &color, IGRAPH_BLISS_F, &info);
    if (strcmp(info.group_size, "1") != 0) {
        printf("Bliss automorphism count with color failed: ring1.\n");
    }
    igraph_free(info.group_size);

    // There's only one automorphism for this coloured graph, so the generating set is empty.
    igraph_automorphism_group_bliss(&ring1, &color, &generators, IGRAPH_BLISS_F, NULL);
    if (igraph_vector_int_list_size(&generators) != 0) {
        printf("Bliss automorphism generators failed with colored graph.\n");
    }

    igraph_vector_int_list_destroy(&generators);

    igraph_vector_int_destroy(&color);

    igraph_vector_int_destroy(&perm);

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);
    igraph_destroy(&directed_ring);
}

void test_bug_995(void) {
    igraph_t g1, g2;
    igraph_bool_t result;

    igraph_small(&g1, 3, 0, 0, 1, 1, 2, 2, 2, -1);
    igraph_small(&g2, 3, 0, 0, 1, 1, 2, 1, 1, -1);

    igraph_isomorphic(&g1, &g2, &result);
    if (result) {
        printf("igraph_isomorphic() failed with loop edges, see bug #995\n");
    }

    igraph_destroy(&g1);
    igraph_destroy(&g2);
}

int main(void) {

    igraph_rng_seed(igraph_rng_default(), 293847); /* make tests deterministic */

    test3();
    test4();
    test_bliss();
    test_bug_995();

    VERIFY_FINALLY_STACK();

    return 0;
}
