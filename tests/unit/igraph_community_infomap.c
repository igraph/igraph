/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"


void gsummary(const igraph_t * g) {
    printf("|V|=%d |E|=%d directed=%d\n", (int)igraph_vcount(g), (int)igraph_ecount(g), (int)igraph_is_directed(g));
}

void show_results(igraph_vector_t * membership, igraph_real_t codelength) {
    int i;
    printf("Codelength: %0.5f (in %d modules)\n", codelength, (int)igraph_vector_max(membership) + 1 );
    printf("Membership: ");
    for (i = 0; i < igraph_vector_size(membership); i++) {
        printf("%li ", (long)VECTOR(*membership)[i] );
    }
    printf("\n");
}

void show_results_lite(igraph_vector_t * membership, igraph_real_t codelength) {
    int i;
    printf("Codelength: %0.5f (in %d modules)\n", codelength, (int)igraph_vector_max(membership) + 1 );
    printf("Membership (1/100 of vertices): ");
    for (i = 0; i < igraph_vector_size(membership); i += 100) {
        printf("%li ", (long)VECTOR(*membership)[i] );
    }
    printf("\n");
}

igraph_real_t infomap_weighted_test(const igraph_t * g, const igraph_vector_t *weights, igraph_bool_t smoke_test) {
    igraph_vector_t membership;
    igraph_real_t codelength = 1000;
    igraph_vector_init(&membership, 0);

    igraph_community_infomap(/*in */ g, /*e_weight=*/ weights, NULL, /*nb_trials=*/5,
                                     /*out*/ &membership, &codelength);
    if (!smoke_test) {
        if (igraph_vcount(g) > 500) {
            show_results_lite(&membership, codelength);
        } else {
            show_results(&membership, codelength);
        }
    }

    igraph_vector_destroy(&membership);

    return codelength;
}


igraph_real_t infomap_test(const igraph_t * g, igraph_bool_t smoke_test) {
    return infomap_weighted_test(g, 0, smoke_test);
}


int main() {
    igraph_t g;
    igraph_vector_t weights;
    igraph_real_t codelength;
    FILE *wikt;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Two triangles connected by one edge */
    printf("# Two triangles connected by one edge\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 0,
                 3, 4, 4, 5, 5, 3,
                 0, 5,
                 -1);
    infomap_test(&g, /* smoke_test = */ 0);
    igraph_destroy(&g);

    /* Two 4-cliques with one commun vertex (vertex 3) */
    printf("# Two 4-cliques (0123 and 4567) connected by two edges (0-4 and 1-5)\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1,  0, 2,  0, 3,  1, 2,  1, 3,  2, 3, /* 4-clique 0,1,2,3 */
                 7, 4,  7, 5,  7, 6,  4, 5,  4, 6,  5, 6, /* 4-clique 4,5,6,7 */
                 0, 4,  1, 5, /* 8, 0, 8, 4, */
                 -1);
    infomap_test(&g, /* smoke_test = */ 0);

    printf("# Two 4-cliques (0123 and 4567) connected by two edges (0-4 and 1-5)\n");
    igraph_add_edge(&g, 0, 4);
    igraph_add_edge(&g, 1, 5);
    infomap_test(&g, /* smoke_test = */ 0);
    igraph_destroy(&g);

    /* Zachary Karate club -- this is just a quick smoke test */
    printf("# Zachary Karate club\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0,  1,  0,  2,  0,  3,  0,  4,  0,  5, /* 0,  5, 0,  5, 0,  5, */
                 0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
                 0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
                 0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
                 1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
                 2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
                 2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
                 4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
                 6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
                 13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
                 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
                 22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
                 23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
                 25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
                 28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
                 31, 32, 31, 33, 32, 33,
                 -1);
    infomap_test(&g, /* smoke_test = */ 0);
    igraph_destroy(&g);

    /* Flow.net that come in infomap_dir.tgz  */
    printf("# Flow (from infomap_dir.tgz)\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1,     1, 2,    2, 3,    3, 0,    1, 4,
                 4, 5,     5, 6,    6, 7,    7, 4,    5, 8,
                 8, 9,     9, 10,  10, 11,  11, 8,    9, 12,
                 12, 13,  13, 14,  14, 15,  15, 12,  13, 0,
                 -1);
    infomap_test(&g, /* smoke_test = */ 0);
    igraph_destroy(&g);

    /* MultiphysChemBioEco40W_weighted_dir.net */
    printf("# MultiphysChemBioEco40W_weighted_dir.net (from infomap_dir.tgz)\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 1, 0,  2, 0,  3, 0,  4, 0,  5, 0,  6, 0,  7, 0,
                 8, 0,  9, 0,  16, 0,  18, 0,  0, 1,  2, 1,  3, 1,
                 5, 1,  6, 1,  7, 1,  9, 1,  10, 1,  16, 1,  18, 1,
                 0, 2,  3, 2,  4, 2,  5, 2,  6, 2,  7, 2,  0, 3,
                 1, 3,  2, 3,  4, 3,  5, 3,  6, 3,  7, 3,  8, 3,
                 9, 3,  10, 3,  11, 3,  13, 3,  14, 3,  16, 3,  17, 3,
                 18, 3,  19, 3,  26, 3,  30, 3,  1, 4,  3, 4,  5, 4,
                 6, 4,  13, 4,  18, 4,  0, 5,  1, 5,  2, 5,  3, 5,
                 6, 5,  7, 5,  9, 5,  1, 6,  3, 6,  7, 6,  9, 6,
                 16, 6,  0, 7,  1, 7,  2, 7,  3, 7,  5, 7,  6, 7,
                 9, 7,  3, 8,  5, 8,  3, 9,  7, 9,  12, 10,  13, 10,
                 14, 10,  15, 10,  16, 10,  17, 10,  18, 10,  19, 10,
                 21, 10,  3, 11,  18, 11,  10, 12,  14, 12,  16, 12,
                 17, 12,  18, 12,  3, 13,  10, 13,  14, 13,  16, 13,
                 10, 14,  12, 14,  13, 14,  15, 14,  16, 14,  17, 14,
                 18, 14,  10, 15,  14, 15,  18, 15,  0, 16,  2, 16,
                 3, 16,  6, 16,  10, 16,  12, 16,  13, 16,  14, 16,
                 17, 16,  18, 16,  10, 17,  12, 17,  14, 17,  18, 17,
                 3, 18,  10, 18,  12, 18,  14, 18,  15, 18,  16, 18,
                 17, 18,  19, 18,  21, 18,  11, 19,  16, 19,  17, 19,
                 16, 20,  18, 20,  21, 20,  22, 20,  23, 20,  24, 20,
                 25, 20,  26, 20,  27, 20,  28, 20,  29, 20,  3, 21,
                 14, 21,  18, 21,  20, 21,  22, 21,  23, 21,  24, 21,
                 25, 21,  26, 21,  27, 21,  28, 21,  29, 21,  35, 21,
                 36, 21,  38, 21,  18, 22,  20, 22,  21, 22,  23, 22,
                 24, 22,  25, 22,  26, 22,  27, 22,  29, 22,  3, 23,
                 20, 23,  21, 23,  22, 23,  24, 23,  25, 23,  26, 23,
                 27, 23,  28, 23,  29, 23,  35, 23,  38, 23,  39, 23,
                 20, 24,  21, 24,  23, 24,  25, 24,  26, 24,  27, 24,
                 28, 24,  29, 24,  9, 25,  20, 25,  21, 25,  22, 25,
                 23, 25,  24, 25,  26, 25,  27, 25,  28, 25,  29, 25,
                 18, 26,  20, 26,  21, 26,  22, 26,  23, 26,  25, 26,
                 27, 26,  28, 26,  29, 26,  30, 26,  32, 26,  35, 26,
                 36, 26,  38, 26,  39, 26,  3, 27,  14, 27,  20, 27,
                 21, 27,  22, 27,  23, 27,  24, 27,  25, 27,  26, 27,
                 28, 27,  29, 27,  38, 27,  3, 28,  18, 28,  20, 28,
                 21, 28,  23, 28,  24, 28,  25, 28,  26, 28,  27, 28,
                 29, 28,  35, 28,  14, 29,  16, 29,  18, 29,  20, 29,
                 21, 29,  22, 29,  23, 29,  24, 29,  25, 29,  26, 29,
                 27, 29,  28, 29,  31, 30,  32, 30,  33, 30,  34, 30,
                 35, 30,  36, 30,  38, 30,  39, 30,  30, 31,  32, 31,
                 34, 31,  36, 31,  30, 32,  34, 32,  35, 32,  36, 32,
                 30, 33,  32, 33,  34, 33,  35, 33,  36, 33,  38, 33,
                 30, 34,  31, 34,  32, 34,  33, 34,  35, 34,  36, 34,
                 38, 34,  39, 34,  26, 35,  30, 35,  32, 35,  33, 35,
                 34, 35,  36, 35,  38, 35,  39, 35,  30, 36,  34, 36,
                 35, 36,  38, 36,  39, 36,  34, 37,  26, 38,  30, 38,
                 32, 38,  33, 38,  34, 38,  35, 38,  36, 38,  39, 38,
                 26, 39,  30, 39,  33, 39,  34, 39,  35, 39,  36, 39,
                 38, 39,
                 -1);
    igraph_vector_init_real(&weights, 306,
                            5.0,  3.0,  130.0,  4.0,  15.0,  9.0,
                            7.0,  1.0,  1.0,  3.0,  1.0,  1.0,
                            1.0,  34.0,  38.0,  2.0,  23.0,  1.0,
                            1.0,  3.0,  2.0,  2.0,  16.0,  1.0,
                            3.0,  1.0,  3.0,  63.0,  92.0,  72.0,
                            25.0,  447.0,  121.0,  65.0,  4.0,  16.0,
                            35.0,  1.0,  19.0,  1.0,  78.0,  1.0,
                            45.0,  1.0,  3.0,  1.0,  1.0,  25.0,
                            1.0,  3.0,  1.0,  1.0,  3.0,  36.0,
                            19.0,  136.0,  41.0,  96.0,  1.0,  7.0,
                            26.0,  1.0,  2.0,  2.0,  3.0,  2.0,  2.0,
                            23.0,  52.0,  4.0,  1.0,  2.0,  1.0,  3.0,
                            1.0,  11.0,  2.0,  17.0,  1.0,  5.0,  18.0,
                            86.0,  5.0,  1.0,  1.0,  1.0,  6.0,  1.0,
                            2.0,  2.0,  20.0,  4.0,  5.0,  1.0,  5.0,
                            12.0,  4.0,  1.0,  1.0,  4.0,  9.0,  40.0,
                            2.0,  1.0,  4.0,  1.0,  1.0,  48.0,  2.0,
                            18.0,  1.0,  7.0,  2.0,  2.0,  53.0,  25.0,
                            9.0,  1.0,  23.0,  8.0,  62.0,  29.0,  35.0,
                            4.0,  34.0,  35.0,  3.0,  1.0,  24.0,  1.0,
                            6.0,  2.0,  2.0,  22.0,  7.0,  2.0,  5.0,
                            14.0,  3.0,  28.0,  14.0,  20.0,  3.0,  1.0,
                            5.0,  77.0,  20.0,  25.0,  35.0,  55.0,  35.0,
                            115.0,  68.0,  105.0,  2.0,  2.0,  2.0,  4.0,
                            2.0,  17.0,  12.0,  3.0,  3.0,  11.0,  10.0,
                            7.0,  2.0,  12.0,  31.0,  11.0,  5.0,  11.0,
                            65.0,  39.0,  17.0,  26.0,  3.0,  4.0,  2.0,
                            3.0,  6.0,  4.0,  8.0,  1.0,  7.0,  7.0,
                            6.0,  1.0,  39.0,  42.0,  9.0,  6.0,  9.0,
                            5.0,  45.0,  43.0,  26.0,  1.0,  2.0,  6.0,
                            2.0,  15.0,  3.0,  9.0,  2.0,  1.0,  1.0,
                            1.0,  4.0,  2.0,  9.0,  2.0,  1.0,  2.0,
                            28.0,  80.0,  10.0,  18.0,  13.0,  17.0,
                            28.0,  40.0,  76.0,  1.0,  2.0,  1.0,  11.0,
                            37.0,  5.0,  11.0,  14.0,  4.0,  14.0,  10.0,
                            1.0,  1.0,  1.0,  1.0,  41.0,  121.0,  6.0,
                            21.0,  12.0,  30.0,  6.0,  141.0,  43.0,  2.0,
                            12.0,  6.0,  35.0,  10.0,  7.0,  2.0,  12.0,
                            6.0,  2.0,  11.0,  1.0,  7.0,  6.0,  5.0,  3.0,
                            1.0,  2.0,  1.0,  1.0,  1.0,  1.0,  67.0,  9.0,
                            9.0,  11.0,  10.0,  21.0,  7.0,  12.0,  9.0,
                            16.0,  7.0,  4.0,  11.0,  17.0,  37.0,  32.0,
                            9.0,  2.0,  2.0,  5.0,  4.0,  2.0,  7.0,  3.0,
                            3.0,  5.0,  8.0,  14.0,  3.0,  38.0,  3.0,  9.0,
                            2.0,  8.0,  21.0,  18.0,  58.0);
    infomap_weighted_test(&g, &weights, /* smoke_test = */ 0);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Wiktionary English verbs -- this one is a bit flaky, igraph reports
     * 1444 or 1445 modules, depending on the platform, so let's just run it as
     * a quick smoke test but don't check the results too thoroughly as some
     * changes are expected. We only check the codelength of the partition,
     * this is more reliable. */
    printf("# Wiktionary english verbs (synonymy 2008)\n");
    wikt = fopen("wikti_en_V_syn.elist", "r");
    igraph_read_graph_edgelist(&g, wikt, 0, 0);
    fclose(wikt);
    gsummary(&g);
    codelength = infomap_test(&g, /* smoke_test = */ 1);
    if (fabs(codelength - 5.708) >= 1e-3) {
        printf("Codelength was %0.5f, expected %0.5f\n", codelength, 5.708);
    } else {
        printf("Codelength OK.\n");
    }
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
