/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

void show_results(igraph_t *g, igraph_vector_t *mod, igraph_matrix_t *merges,
                  igraph_vector_t *membership, FILE* f) {
    long int i = 0;
    igraph_vector_t our_membership;

    igraph_vector_init(&our_membership, 0);

    if (mod != 0) {
        i = igraph_vector_which_max(mod);
        fprintf(f, "Modularity:  %f\n", VECTOR(*mod)[i]);
    } else {
        fprintf(f, "Modularity:  ---\n");
    }

    if (membership != 0) {
        igraph_vector_update(&our_membership, membership);
    } else if (merges != 0) {
        igraph_community_to_membership(merges, igraph_vcount(g), i, &our_membership, 0);
    }

    printf("Membership: ");
    for (i = 0; i < igraph_vector_size(&our_membership); i++) {
        printf("%li ", (long int)VECTOR(our_membership)[i]);
    }
    printf("\n");

    igraph_vector_destroy(&our_membership);
}

int main() {
    igraph_t g;
    igraph_vector_t modularity, weights, membership;
    igraph_matrix_t merges;

    igraph_vector_init(&modularity, 0);
    igraph_matrix_init(&merges, 0, 0);
    igraph_vector_init(&weights, 0);
    igraph_vector_init(&membership, 0);

    /* Simple unweighted graph */
    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);

    /* Same simple graph, with uniform edge weights */
    igraph_vector_resize(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 2);
    igraph_community_fastgreedy(&g, &weights, &merges, &modularity,
                                /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Simple nonuniform weighted graph, with and without weights */
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 2, 4, 2, 5, 3, 4, 3, 5, 4, 5, -1);
    igraph_vector_resize(&weights, 8);
    igraph_vector_fill(&weights, 1);
    VECTOR(weights)[0] = 10;
    VECTOR(weights)[1] = 10;
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_community_fastgreedy(&g, &weights, &merges, &modularity,
                                /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Zachary Karate club */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
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
    igraph_community_fastgreedy(&g, 0, &merges, &modularity,
                                /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Simple disconnected graph with isolates */
    igraph_small(&g, 9, IGRAPH_UNDIRECTED,
                 0,  1,  0,  2,  0,  3,  1,  2,  1,  3,  2,  3,
                 4,  5,  4,  6,  4,  7,  5,  6,  5,  7,  6,  7,
                 -1);
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Disjoint union of two rings */
    igraph_small(&g, 20, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 0, 9,
                 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 10, 19, -1);
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Completely empty graph */
    igraph_small(&g, 10, IGRAPH_UNDIRECTED, -1);
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Ring graph with loop edges */
    igraph_small(&g, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 2, 2, -1);
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Regression test -- graph with two vertices and two edges */
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 0, 1, 1, -1);
    igraph_community_fastgreedy(&g, 0, &merges, &modularity, /*membership=*/ 0);
    show_results(&g, &modularity, &merges, 0, stdout);
    igraph_destroy(&g);

    /* Regression test -- asking for optimal membership vector but not
     * providing a modularity vector */
    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);
    igraph_community_fastgreedy(&g, 0, &merges, 0, &membership);
    show_results(&g, 0, &merges, &membership, stdout);
    igraph_destroy(&g);

    /* Regression test -- asking for optimal membership vector but not
     * providing a merge matrix */
    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);
    igraph_community_fastgreedy(&g, 0, 0, &modularity, &membership);
    show_results(&g, &modularity, 0, &membership, stdout);

    /* Regression test -- asking for optimal membership vector but not
     * providing a merge matrix or a modularity vector */
    igraph_community_fastgreedy(&g, 0, 0, 0, &membership);
    show_results(&g, 0, 0, &membership, stdout);
    igraph_destroy(&g);

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&modularity);
    igraph_vector_destroy(&weights);
    igraph_matrix_destroy(&merges);

    return 0;
}
