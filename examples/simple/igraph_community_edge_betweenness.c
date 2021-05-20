/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int igraph_vector_between(const igraph_vector_t* v, const igraph_vector_t* lo,
                          const igraph_vector_t* hi) {
    return igraph_vector_all_le(lo, v) && igraph_vector_all_ge(hi, v);
}

void test_unweighted() {
    igraph_t g;
    igraph_vector_t edges, eb;
    long int i;
    long int no_of_edges;

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

    igraph_vector_init(&edges, 0);
    igraph_vector_init(&eb, 0);
    igraph_community_edge_betweenness(&g, &edges, &eb, 0 /*merges */,
                                      0 /*bridges */, /*modularity=*/ 0,
                                      /*membership=*/ 0,
                                      IGRAPH_UNDIRECTED,
                                      /*weights=*/ 0);

    no_of_edges = igraph_ecount(&g);
    for (i = 0; i < no_of_edges; i++) {
        printf("%li ", (long int)VECTOR(edges)[i]);
    }
    printf("\n");

    for (i = 0; i < no_of_edges; i++) {
        printf("%.2f ", VECTOR(eb)[i]);
    }
    printf("\n");

    /* Try it once again without storage space for edges */
    igraph_community_edge_betweenness(&g, 0, &eb, 0 /*merges */,
                                      0 /*bridges */, /*modularity=*/ 0,
                                      /*membership=*/ 0,
                                      IGRAPH_UNDIRECTED,
                                      /*weights=*/ 0);
    for (i = 0; i < no_of_edges; i++) {
        printf("%.2f ", VECTOR(eb)[i]);
    }
    printf("\n");

    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&edges);
    igraph_destroy(&g);
}

#define EPS 1e-4

void test_weighted() {
    igraph_t g;
    igraph_vector_t edges, eb, weights;
    igraph_real_t weights_array[] = { 4, 1, 3, 2, 5, 8, 6, 7 };

    igraph_real_t edges_array1[] = { 2, 3, 0, 1, 4, 7, 5, 6 };
    igraph_real_t edges_array2[] = { 2, 3, 6, 5, 0, 1, 4, 7 };
    igraph_real_t eb_array1_lo[] = { 4, 5, 3 + 1 / 3.0 - EPS, 4, 2.5, 4, 1, 1 };
    igraph_real_t eb_array1_hi[] = { 4, 5, 3 + 1 / 3.0 + EPS, 4, 2.5, 4, 1, 1 };
    igraph_real_t eb_array2_lo[] = { 4, 5, 3 + 1 / 3.0 - EPS, 6, 1.5, 2, 1, 1 };
    igraph_real_t eb_array2_hi[] = { 4, 5, 3 + 1 / 3.0 + EPS, 6, 1.5, 2, 1, 1 };

    igraph_vector_t edges_sol1, edges_sol2, eb_sol1_lo, eb_sol1_hi, eb_sol2_lo, eb_sol2_hi;

    igraph_vector_view(&edges_sol1, edges_array1,
                       sizeof(edges_array1) / sizeof(double));
    igraph_vector_view(&edges_sol2, edges_array2,
                       sizeof(edges_array2) / sizeof(double));
    igraph_vector_view(&eb_sol1_lo, eb_array1_lo, sizeof(eb_array1_lo) / sizeof(double));
    igraph_vector_view(&eb_sol2_lo, eb_array2_lo, sizeof(eb_array2_lo) / sizeof(double));
    igraph_vector_view(&eb_sol1_hi, eb_array1_hi, sizeof(eb_array1_hi) / sizeof(double));
    igraph_vector_view(&eb_sol2_hi, eb_array2_hi, sizeof(eb_array2_hi) / sizeof(double));

    /* Small graph as follows: A--B--C--A, A--D--E--A, B--D, C--E */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 2, 4, 3, 4, -1);
    igraph_vector_view(&weights, weights_array, igraph_ecount(&g));

    igraph_vector_init(&edges, 0);
    igraph_vector_init(&eb, 0);
    igraph_community_edge_betweenness(&g, &edges, &eb, 0 /*merges */,
                                      0 /*bridges */, /*modularity=*/ 0,
                                      /*membership=*/ 0,
                                      IGRAPH_UNDIRECTED,
                                      &weights);

    if (!igraph_vector_all_e(&edges_sol1, &edges) &&
        !igraph_vector_all_e(&edges_sol2, &edges)) {
        printf("Error, edges vector was: \n");
        igraph_vector_print(&edges);
        exit(2);
    }
    if (!igraph_vector_between(&eb, &eb_sol1_lo, &eb_sol1_hi) &&
        !igraph_vector_between(&eb, &eb_sol2_lo, &eb_sol2_hi)) {
        printf("Error, eb vector was: \n");
        igraph_vector_print(&eb);
        exit(2);
    }

    /* Try it once again without storage space for edges */
    igraph_community_edge_betweenness(&g, 0, &eb, 0 /*merges */,
                                      0 /*bridges */, /*modularity=*/ 0,
                                      /*membership=*/ 0,
                                      IGRAPH_UNDIRECTED,
                                      &weights);

    if (!igraph_vector_between(&eb, &eb_sol1_lo, &eb_sol1_hi) &&
        !igraph_vector_between(&eb, &eb_sol2_lo, &eb_sol2_hi)) {
        printf("Error, eb vector was: \n");
        igraph_vector_print(&eb);
        exit(2);
    }

    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&edges);
    igraph_destroy(&g);
}

void test_zero_edge_graph() {
    igraph_t g;
    igraph_vector_t eb;
    igraph_vector_t res;

    igraph_full(&g, 1, 0, 0);
    igraph_vector_init(&res, igraph_ecount(&g));
    igraph_vector_init(&eb, igraph_ecount(&g));

    igraph_community_edge_betweenness(&g,
        &res, // result
        &eb, // edge_betweenness result
        NULL, // merges result
        NULL, // bridges
        NULL, // modularity
        NULL, // membership
        IGRAPH_UNDIRECTED, // directed
        NULL // weights
        );

    igraph_vector_destroy(&eb);
    printf("No crash\n");
    igraph_vector_destroy(&res);
    igraph_destroy(&g);
}

int main() {
    test_unweighted();
    test_weighted();
    test_zero_edge_graph();
    return 0;
}
