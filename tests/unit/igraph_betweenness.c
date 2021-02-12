/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main() {

    igraph_t g;
    igraph_vector_t bet, bet2, weights, edges;

    igraph_real_t cutoff = 0.0;

    igraph_real_t nontriv[] = { 0, 19, 0, 16, 0, 20, 1, 19, 2, 5, 3, 7, 3, 8,
                                4, 15, 4, 11, 5, 8, 5, 19, 6, 7, 6, 10, 6, 8,
                                6, 9, 7, 20, 9, 10, 9, 20, 10, 19,
                                11, 12, 11, 20, 12, 15, 13, 15,
                                14, 18, 14, 16, 14, 17, 15, 16, 17, 18
                              };

    igraph_real_t nontriv_weights[] = { 0.5249, 1, 0.1934, 0.6274, 0.5249,
                                        0.0029, 0.3831, 0.05, 0.6274, 0.3831,
                                        0.5249, 0.0587, 0.0579, 0.0562, 0.0562,
                                        0.1934, 0.6274, 0.6274, 0.6274, 0.0418,
                                        0.6274, 0.3511, 0.3511, 0.1486, 1, 1,
                                        0.0711, 0.2409
                                      };

    igraph_real_t nontriv_res[] = { 20, 0, 0, 0, 0, 19, 80, 85, 32, 0, 10,
                                    75, 70, 0, 36, 81, 60, 0, 19, 19, 86
                                  };

    /*******************************************************/

    printf("BA graph\n");
    printf("==========================================================\n");
    igraph_barabasi_game(/* graph= */    &g,
                                         /* n= */        1000,
                                         /* power= */    1,
                                         /* m= */        3,
                                         /* outseq= */   0,
                                         /* outpref= */  0,
                                         /* A= */        1,
                                         /* directed= */ 0,
                                         /* algo= */     IGRAPH_BARABASI_BAG,
                                         /* start_from= */ 0);

    igraph_simplify(&g, /* multiple= */ 1, /* loops= */ 1, /*edge_comb=*/ 0);

    igraph_vector_init(&bet, 0);

    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ 0,
            /* cutoff=    */ 2);

    igraph_vector_destroy(&bet);
    igraph_destroy(&g);

    printf("\nTree\n");
    printf("==========================================================\n");
    igraph_tree(&g, 20000, 10, IGRAPH_TREE_UNDIRECTED);

    igraph_vector_init(&bet, 0);

    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ 0,
            /* cutoff=    */ 3);

    printf("Max betweenness: %f\n", igraph_vector_max(&bet));

    igraph_vector_init(&bet2, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1.0);

    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ &weights,
            /* cutoff=    */ 3);

    IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nNon-trivial weighted graph\n");
    printf("==========================================================\n");
    igraph_vector_view(&edges, nontriv, sizeof(nontriv) / sizeof(igraph_real_t));
    igraph_create(&g, &edges, 0, /* directed= */ 0);
    igraph_vector_view(&weights, nontriv_weights,
                       sizeof(nontriv_weights) / sizeof(igraph_real_t));
    igraph_vector_init(&bet, 0);

    igraph_betweenness(/*graph=*/ &g, /*res=*/ &bet, /*vids=*/ igraph_vss_all(),
                                  /*directed=*/0, /*weights=*/ &weights);

    printf("Max betweenness: %f\n", igraph_vector_max(&bet));

    igraph_vector_view(&bet2, nontriv_res,
                       sizeof(nontriv_res) / sizeof(igraph_real_t));

    IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));

    igraph_vector_destroy(&bet);
    igraph_destroy(&g);


    printf("\nCorner case cutoff 0.0\n");
    printf("==========================================================\n");
    igraph_tree(&g, 20, 3, IGRAPH_TREE_UNDIRECTED);

    /* unweighted */
    igraph_vector_init(&bet, 0);
    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ 0,
            /* cutoff=    */ 0);

    igraph_vector_init(&bet2, 0);
    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ 0,
            /* cutoff=    */ -1);

    igraph_vector_print(&bet);
    igraph_vector_print(&bet2);

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);

    /* weighted */
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 2.0);

    igraph_vector_init(&bet, 0);
    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ &weights,
            /* cutoff=    */ 0);

    igraph_vector_init(&bet2, 0);
    igraph_betweenness_cutoff(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* weights=   */ &weights,
            /* cutoff=    */ -1);

    igraph_vector_print(&bet);
    igraph_vector_print(&bet2);

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nSingle path graph\n");
    printf("==========================================================\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED,
                            0, 1,
                            1, 2,
                            2, 3,
                            3, 4, -1);
    igraph_vector_init(&bet, igraph_vcount(&g));
    igraph_vector_init(&bet2, igraph_vcount(&g));
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1);

    for (cutoff = -1.0; cutoff < 5.0; cutoff += 1)
    {
        printf("Cutoff %.0f\n", cutoff);
        printf("Unweighted\n");
        igraph_betweenness_cutoff(&g, &bet,
                                  igraph_vss_all(), IGRAPH_UNDIRECTED,
                /* weights */ NULL,
                /* cutoff */ cutoff);
        igraph_vector_print(&bet);

        printf("Weighted\n");
        igraph_betweenness_cutoff(&g, &bet2,
                                  igraph_vss_all(), IGRAPH_UNDIRECTED,
                /* weights */ &weights,
                /* cutoff */ cutoff);
        igraph_vector_print(&bet2);
        printf("\n");

        IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));
    }

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nCycle graph\n");
    printf("==========================================================\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED,
                            0, 1,
                            0, 2,
                            1, 3,
                            2, 3, -1);
    igraph_vector_init(&bet, igraph_vcount(&g));
    igraph_vector_init(&bet2, igraph_vcount(&g));
    igraph_vector_init(&weights, igraph_ecount(&g));
    VECTOR(weights)[0] = 1.01;
    VECTOR(weights)[1] = 2;
    VECTOR(weights)[2] = 0.99;
    VECTOR(weights)[3] = 2;

    for (cutoff = -1.0; cutoff < 5.0; cutoff += 1)
    {
        printf("Cutoff %.0f\n", cutoff);
        printf("Unweighted\n");
        igraph_betweenness_cutoff(&g, &bet,
                                  igraph_vss_all(), IGRAPH_UNDIRECTED,
                /* weights */ NULL,
                /* cutoff */ cutoff);
        igraph_vector_print(&bet);

        printf("Weighted\n");
        igraph_betweenness_cutoff(&g, &bet2,
                                  igraph_vss_all(), IGRAPH_UNDIRECTED,
                /* weights */ &weights,
                /* cutoff */ cutoff);
        igraph_vector_print(&bet2);
        printf("\n");
    }

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nNull graph\n");
    printf("==========================================================\n");

    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    igraph_vector_init(&bet, 3); /* purposefully larger than zero, as igraph_betweenness must resize it */
    igraph_betweenness(&g, &bet, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL);
    print_vector(&bet);

    igraph_vector_destroy(&bet);
    igraph_destroy(&g);

    printf("\nEmpty graph\n");
    printf("==========================================================\n");

    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);
    igraph_vector_init(&bet, 0);
    igraph_betweenness(&g, &bet, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL);
    print_vector(&bet);

    igraph_vector_destroy(&bet);
    igraph_destroy(&g);

    printf("\n37x37 grid graph\n");
    printf("==========================================================\n");

    {
        igraph_vector_t dims;

        igraph_vector_init(&dims, 2);
        VECTOR(dims)[0] = 37;
        VECTOR(dims)[1] = 37;

        igraph_lattice(&g, &dims, 1, IGRAPH_UNDIRECTED, 0, 0);

        igraph_vector_init(&bet, 0);
        igraph_betweenness(&g, &bet, igraph_vss_all(), IGRAPH_UNDIRECTED, NULL);
        printf("Max betweenness: %f\n", igraph_vector_max(&bet));

        igraph_vector_destroy(&bet);
        igraph_destroy(&g);
        igraph_vector_destroy(&dims);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}

