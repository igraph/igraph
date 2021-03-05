/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2021  Gabor Csardi <csardi.gabor@gmail.com>
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


/* https://github.com/igraph/igraph/issues/950 */
void test_bug950_edge() {
    /* Testing the case of weighted graphs with multiple alternate
     * paths to the same node with slightly different weights due to
     * floating point inaccuracies. */
    igraph_t g;
    igraph_vector_t eb;
    igraph_vector_t weights;
    igraph_integer_t from, to;
    long int no_of_edges, i;

    igraph_full(&g, 6, 0, 0);
    no_of_edges = igraph_ecount(&g);

    igraph_vector_init(&weights, no_of_edges);

    for (i = 0; i < no_of_edges; i++) {
        igraph_edge(&g, i, &from, &to);
        if((from < 3 && to < 3) || (from >= 3 && to >= 3))
            VECTOR(weights)[i] = 1;
        else
            VECTOR(weights)[i] = 0.1;
    }

    printf("\n testing bug 950 source and targets = vss_all() \n");
    printf("==========================================================\n");
    igraph_vector_init(&eb, 0);
    igraph_edge_betweenness_subset(&g, &eb, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, igraph_vss_all(), igraph_vss_all(), &weights);
    print_vector(&eb);

    igraph_vector_destroy(&eb);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
}

int main() {
    igraph_t g;
    igraph_es_t es;
    igraph_vector_t eb, node_vec, source_vec, target_vec, bet, bet2, weights, edges;
    igraph_vs_t vs, vs_source, vs_target;

    /* edge betweenness test */

    {
        /* We use igraph_create() instead of igraph_small() as some MSVC versions
           will choke on an overlong argument list with "internal error C1001". */
        igraph_real_t edge_array[] = {
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
            31, 32, 31, 33, 32, 33
        };

        igraph_create(&g, igraph_vector_view(&edges, edge_array, sizeof(edge_array) / sizeof(igraph_real_t)), 0, IGRAPH_UNDIRECTED);
        igraph_vector_init_seq(&source_vec, 0, 32);
        igraph_vs_vector(&vs_source, &source_vec);
        igraph_vector_init_seq(&target_vec, 1, 33);
        igraph_vs_vector(&vs_target, &target_vec);
        igraph_vector_init(&eb, 0);
        printf("==========================================================\n");
        printf("\nlarge unweighted graph edge betweenness\n");
        igraph_edge_betweenness_subset(&g, &eb, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, vs_source, vs_target, /*weights=*/ 0);
        print_vector(&eb);
        igraph_vector_destroy(&eb);
        igraph_destroy(&g);
        igraph_vs_destroy(&vs_source);
        igraph_vector_destroy(&source_vec);
        igraph_vs_destroy(&vs_target);
        igraph_vector_destroy(&target_vec);
    }

    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 1, 4, -1);

    igraph_vector_init(&eb, 0);
    igraph_vector_init_seq(&node_vec, 0, 3);
    igraph_es_vector(&es, &node_vec);
    igraph_vector_init_seq(&target_vec, 1, 4);
    igraph_vs_vector(&vs_target, &target_vec);
    igraph_edge_betweenness_subset(&g, &eb, es, IGRAPH_UNDIRECTED, igraph_vss_all(), vs_target, /*weights=*/ 0);
    printf("==========================================================\n");
    printf("\n small unweighted edge betweenness\n");
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_es_destroy(&es);
    igraph_vector_destroy(&node_vec);
    igraph_vs_destroy(&vs_target);
    igraph_vector_destroy(&target_vec);
    igraph_destroy(&g);


    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 3, 1, 2, 1, 4, 2, 5, 3, 4, 3, 6, 4, 5, 4, 7, 5, 8,
                 6, 7, 7, 8, -1);
    igraph_vector_init(&eb, 0);
    igraph_vector_init_seq(&source_vec, 1, 8);
    igraph_vs_vector(&vs_source, &source_vec);
    igraph_edge_betweenness_subset(&g, &eb, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, vs_source, igraph_vss_all(), /*weights=*/ 0);
    printf("==========================================================\n");
    printf("\n unweighted small graph edge betweenness with subset of sources");
    print_vector(&eb);
    igraph_vector_destroy(&eb);
    igraph_vs_destroy(&vs_source);
    igraph_vector_destroy(&source_vec);
    igraph_destroy(&g);

    test_bug950_edge();

    int x = IGRAPH_FINALLY_STACK_SIZE();


    /*******************************************************/

    /* vertex betweenness test */
    
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
    igraph_vs_seq(&vs_source, 0, 500);
    igraph_vs_seq(&vs_target, 500, 999);

    igraph_betweenness_subset(/* graph=     */ &g,
        /* res=       */ &bet,
        /* vids=      */ igraph_vss_all(),
        /* directed = */ 0,
        /* sources = */ vs_source,
        /* target = */ vs_target,
        /* weights=   */ 0);

        igraph_vector_destroy(&bet);
        igraph_vs_destroy(&vs_source);
        igraph_vs_destroy(&vs_target);
        igraph_destroy(&g);

    printf("\nTree\n");
    printf("==========================================================\n");
    igraph_tree(&g, 20000, 10, IGRAPH_TREE_UNDIRECTED);

    igraph_vs_seq(&vs_source, 19800, 19999); 
    igraph_vs_seq(&vs_target, 19800, 19999);    
    igraph_vector_init(&bet, 0);

    igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ vs_source,
            /* target = */ vs_target,
            /* weights=   */ 0);

    printf("Max betweenness: %f\n", igraph_vector_max(&bet));

    igraph_vector_init(&bet2, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1.0);

    igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ vs_source,
            /* target = */ vs_target,
            /* weights=   */ &weights);

    IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));


    igraph_vector_destroy(&weights);
    igraph_vs_destroy(&vs_source);
    igraph_vs_destroy(&vs_target);
    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_destroy(&g);

    printf("\nnon-trivial weighted graph using subset algorithm\n");
    printf("==========================================================\n");
    igraph_vector_view(&edges, nontriv, sizeof(nontriv) / sizeof(igraph_real_t));
    igraph_create(&g, &edges, 0, /* directed= */ 0);
    igraph_vector_view(&weights, nontriv_weights,
                       sizeof(nontriv_weights) / sizeof(igraph_real_t));

    igraph_vector_init(&bet, 0);
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ igraph_vss_all(),
            /* target = */ igraph_vss_all(),
            /* weights=   */ &weights);
        igraph_vector_view(&bet2, nontriv_res,
                       sizeof(nontriv_res) / sizeof(igraph_real_t));
    IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));

    igraph_vector_destroy(&bet);
    igraph_destroy(&g);

    printf("\nSingle path graph of subset\n");
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

    for (int i = 0; i<5; i++)
    {
        igraph_vector_init_seq(&node_vec, 0, 4);
        igraph_vector_remove(&node_vec, (long int) i);
        igraph_vs_vector(&vs, &node_vec);
        igraph_vector_init_seq(&source_vec, 0, 4);
        igraph_vector_remove(&source_vec, (long int) i);
        igraph_vs_vector(&vs_source, &source_vec);
        printf("subset without %d\n", i);
        printf("Unweighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ vs,
            /* directed = */ 0,
            /* sources = */ vs_source,
            /* target = */ igraph_vss_all(),
            /* weights=   */ 0);
        igraph_vector_print(&bet);

        printf("Weighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ vs,
            /* directed = */ 0,
            /* sources = */ vs_source,
            /* target = */ igraph_vss_all(),
            /* weights */ &weights);
        igraph_vector_print(&bet2);
        printf("\n");

        IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));
        igraph_vs_destroy(&vs);
        igraph_vector_destroy(&node_vec);
        igraph_vs_destroy(&vs_source);
        igraph_vector_destroy(&source_vec);
    }

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nCycle graph subset\n");
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

    for (int i = 0; i<3; i++)
    {
        igraph_vector_init_seq(&target_vec, 0, 3);
        igraph_vector_remove(&target_vec, (long int) i);
        igraph_vs_vector(&vs_target, &target_vec);
        printf("subset without %d\n", i);
        printf("Unweighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ igraph_vss_all(),
            /* target = */ vs_target,
            /* weights */ NULL);
        igraph_vector_print(&bet);

        printf("Weighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ igraph_vss_all(),
            /* target = */ vs_target,
            /* weights */ &weights);
        igraph_vector_print(&bet2);
        printf("\n");

        igraph_vs_destroy(&vs_target);
        igraph_vector_destroy(&target_vec);
    }

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nEmpty graph\n");
    printf("==========================================================\n");
    
    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);
    igraph_vector_init(&bet, 0);
    igraph_betweenness_subset(&g, &bet, igraph_vss_all(), IGRAPH_UNDIRECTED, igraph_vss_all(), igraph_vss_all(), NULL);
    
    print_vector(&bet);

    igraph_vector_destroy(&bet);
    
    igraph_vector_init(&bet, 0);
    igraph_edge_betweenness_subset(&g, &bet, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, igraph_vss_all(), igraph_vss_all(), NULL);

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
        igraph_vector_init (&bet2, 0);
        igraph_vector_init_seq(&target_vec, 0, (int) igraph_vcount(&g) - 1);
        igraph_vector_remove(&target_vec, (long int)0);
        igraph_vs_vector(&vs_target, &target_vec);
        igraph_vector_init_seq(&source_vec, 0, (int) igraph_vcount(&g) - 1);
        igraph_vector_remove(&source_vec, (long int) 0);
        igraph_vs_vector(&vs_source, &source_vec);
        igraph_betweenness_subset(&g, &bet, igraph_vss_all(), IGRAPH_UNDIRECTED, vs_source, vs_target, NULL);
        printf("Max betweenness: %f\n", igraph_vector_max(&bet));

        igraph_edge_betweenness_subset(&g, &bet2, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, vs_source, vs_target, NULL);
        printf("Max betweenness: %f\n", igraph_vector_max(&bet2));

        igraph_vector_destroy(&bet);
        igraph_vector_destroy(&bet2);
        igraph_destroy(&g);
        igraph_vector_destroy(&dims);
        igraph_vs_destroy(&vs_target);
        igraph_vector_destroy(&target_vec);
        igraph_vs_destroy(&vs_source);
        igraph_vector_destroy(&source_vec);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
