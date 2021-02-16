
#include <igraph.h>
#include <stdio.h>
#include <string.h>

#include "test_utilities.inc"

int random_permutation(igraph_vector_t *vec) {
    /* We just do size(vec) * 2 swaps */
    long int one, two, tmp, i, n = igraph_vector_size(vec);
    for (i = 0; i < 2 * n; i++) {
        one = (double)rand() / RAND_MAX * n;
        two = (double)rand() / RAND_MAX * n;
        tmp = one;
        one = two;
        two = tmp;
    }
    return 0;
}


void test3() {
    int i, j;
    igraph_vector_ptr_t graphs3;

    // Verify that no two 3-vertex graphs of distinct isoclasses are considered isomorphic by Bliss or VF2.

    igraph_vector_ptr_init(&graphs3, 0);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&graphs3, igraph_destroy);

    for (i = 0; i < 16; i++) {
        igraph_t *g;
        g = (igraph_t *) malloc(sizeof(igraph_t));
        igraph_vector_ptr_push_back(&graphs3, g);
        igraph_isoclass_create(g, 3, i, /* directed = */ 1);
    }

    for (i = 0; i < 16; i++)
        for (j = i + 1; j < 16; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_bliss(
                (igraph_t *) VECTOR(graphs3)[i], (igraph_t *) VECTOR(graphs3)[j],
                NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
            if (iso) {
                printf("Bliss failure, 3 vertex directed graphs of isoclass %d and %d are not isomorphic. Bliss reports otherwise.\n", i, j);
            }
        }

    for (i = 0; i < 16; i++)
        for (j = i + 1; j < 16; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_vf2(
                (igraph_t *) VECTOR(graphs3)[i], (igraph_t *) VECTOR(graphs3)[j],
                NULL, NULL, NULL, NULL, &iso, NULL, NULL, NULL, NULL, NULL);
            if (iso) {
                printf("VF2 failure, 3 vertex directed graphs of isoclass %d and %d are not isomorphic. VF2 reports otherwise.\n", i, j);
            }
        }

    igraph_vector_ptr_destroy_all(&graphs3);
}


void test4() {
    int i, j;
    igraph_vector_ptr_t graphs4;

    // Verify that no two 4-vertex graphs of distinct isoclasses are considered isomorphic by Bliss or VF2.

    igraph_vector_ptr_init(&graphs4, 0);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&graphs4, igraph_destroy);

    for (i = 0; i < 218; i++) {
        igraph_t *g;
        g = (igraph_t *) malloc(sizeof(igraph_t));
        igraph_vector_ptr_push_back(&graphs4, g);
        igraph_isoclass_create(g, 4, i, /* directed = */ 1);
    }

    for (i = 0; i < 218; i++)
        for (j = i + 1; j < 218; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_bliss(
                (igraph_t *) VECTOR(graphs4)[i], (igraph_t *) VECTOR(graphs4)[j],
                NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
            if (iso) {
                printf("Bliss failure, 4 vertex directed graphs of isoclass %d and %d are not isomorphic. Bliss reports otherwise.\n", i, j);
            }
        }

    for (i = 0; i < 218; i++)
        for (j = i + 1; j < 218; j++) {
            igraph_bool_t iso;
            igraph_isomorphic_vf2(
                (igraph_t *) VECTOR(graphs4)[i], (igraph_t *) VECTOR(graphs4)[j],
                NULL, NULL, NULL, NULL, &iso, NULL, NULL, NULL, NULL, NULL);
            if (iso) {
                printf("VF2 failure, 4 vertex directed graphs of isoclass %d and %d are not isomorphic. VF2 reports otherwise.\n", i, j);
            }
        }

    igraph_vector_ptr_destroy_all(&graphs4);
}


void test_bliss() {
    igraph_t ring1, ring2, directed_ring;
    igraph_vector_t perm;
    igraph_bool_t iso;
    igraph_bliss_info_t info;
    igraph_vector_int_t color;
    igraph_vector_ptr_t generators;

    igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/1);
    igraph_vector_init_seq(&perm, 0, igraph_vcount(&ring1) - 1);
    random_permutation(&perm);
    igraph_permute_vertices(&ring1, &ring2, &perm);

    igraph_ring(&directed_ring, 100, /* directed= */ 1, /* mutual = */0, /* circular = */1);

    igraph_vector_ptr_init(&generators, 0);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&generators, igraph_vector_destroy);

    igraph_isomorphic_bliss(&ring1, &ring2, NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
    if (! iso) {
        printf("Bliss failed on ring isomorphism.\n");
    }

    igraph_automorphisms(&ring1, NULL, IGRAPH_BLISS_F, &info);
    if (strcmp(info.group_size, "200") != 0) {
        printf("Biss automorphism count failed: ring1.\n");
    }
    igraph_free(info.group_size);

    igraph_automorphisms(&ring2, NULL, IGRAPH_BLISS_F, &info);
    if (strcmp(info.group_size, "200") != 0) {
        printf("Biss automorphism count failed: ring2.\n");
    }
    igraph_free(info.group_size);

    igraph_automorphisms(&directed_ring, NULL, IGRAPH_BLISS_F, &info);
    if (strcmp(info.group_size, "100") != 0) {
        printf("Biss automorphism count failed: directed_ring.\n");
    }
    igraph_free(info.group_size);

    // The follwing test is included so there is at least one call to igraph_automorphism_group
    // in the test suite. However, the generator set returned may depend on the splitting
    // heursitics as well as on the Bliss version. If the test fails, please verify manually
    // that the generating set is valid. For a undirected cycle graph like ring2, there should
    // be two generators: a cyclic permutation and a reversal of the vertex order.
    igraph_automorphism_group(&ring2, NULL, &generators, IGRAPH_BLISS_F, NULL);
    if (igraph_vector_ptr_size(&generators) != 2)
        printf("Bliss automorphism generators may have failed with ring2. "
               "Please verify the generators manually. "
               "Note that the generator set is not guaranteed to be minimal.\n");
    igraph_vector_ptr_free_all(&generators);

    // For a directed ring, the only generator should be a cyclic permutation.
    igraph_automorphism_group(&directed_ring, NULL, &generators, IGRAPH_BLISS_F, NULL);
    if (igraph_vector_ptr_size(&generators) != 1)
        printf("Bliss automorphism generators may have failed with directed_ring. "
               "Please verify the generators manually. "
               "Note that the generator set is not guaranteed to be minimal.\n");
    igraph_vector_ptr_free_all(&generators);

    igraph_vector_int_init_seq(&color, 0, igraph_vcount(&ring1) - 1);

    igraph_automorphisms(&ring1, &color, IGRAPH_BLISS_F, &info);
    if (strcmp(info.group_size, "1") != 0) {
        printf("Biss automorphism count with color failed: ring1.\n");
    }
    igraph_free(info.group_size);

    // There's only one automorphism for this coloured graph, so the generating set is empty.
    igraph_automorphism_group(&ring1, &color, &generators, IGRAPH_BLISS_F, NULL);
    if (igraph_vector_ptr_size(&generators) != 0) {
        printf("Bliss automorphism generators failed with colored graph.\n");
    }

    igraph_vector_ptr_destroy_all(&generators);

    igraph_vector_int_destroy(&color);

    igraph_vector_destroy(&perm);

    igraph_destroy(&ring1);
    igraph_destroy(&ring2);
    igraph_destroy(&directed_ring);
}

void test_bug_995() {
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

int main() {

    srand(293847); /* rand() is used in random_permutation() */

    test3();
    test4();
    test_bliss();
    test_bug_995();

    VERIFY_FINALLY_STACK();

    return 0;
}
