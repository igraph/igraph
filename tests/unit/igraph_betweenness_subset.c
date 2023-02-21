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
#include "test_utilities.h"


int main(void) {
    igraph_t g;
    igraph_vector_int_t edges;
    igraph_vector_t bet, bet2, weights;
    igraph_vector_int_t node_vec, source_vec, target_vec;
    igraph_vs_t vs, vs_source, vs_target;
    igraph_integer_t i, n;

    igraph_integer_t nontriv[] = { 0, 19, 0, 16, 0, 20, 1, 19, 2, 5, 3, 7, 3, 8,
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

    printf("BA graph (smoke test)\n");
    printf("==========================================================\n");
    igraph_barabasi_game(/* graph= */    &g,
                                         /* n= */        1000,
                                         /* power= */    1,
                                         /* m= */        3,
                                         /* outseq= */   0,
                                         /* outpref= */  0,
                                         /* A= */        1,
                                         /* directed= */ IGRAPH_UNDIRECTED,
                                         /* algo= */     IGRAPH_BARABASI_BAG,
                                         /* start_from= */ 0);

    igraph_simplify(&g, /* multiple= */ true, /* loops= */ true, /*edge_comb=*/ NULL);

    igraph_vector_init(&bet, 0);
    igraph_vs_range(&vs_source, 0, 501);
    igraph_vs_range(&vs_target, 500, 1000);

    igraph_betweenness_subset(/* graph=     */ &g,
        /* res=       */ &bet,
        /* vids=      */ igraph_vss_all(),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ vs_source,
        /* target = */ vs_target,
        /* weights=   */ NULL);

    igraph_vector_destroy(&bet);
    igraph_vs_destroy(&vs_source);
    igraph_vs_destroy(&vs_target);
    igraph_destroy(&g);

    printf("\nTree\n");
    printf("==========================================================\n");
    igraph_kary_tree(&g, 11111, 10, IGRAPH_TREE_UNDIRECTED);

    /* We are including the rightmost 200 vertices from the lowermost layer
     * (layer 5) of the tree. These have 20 parents in layer 4, 2 grandparents
     * in layer 3, and a single grand-grandparent in layer 2. None of the
     * shortest paths we consider should therefore pass through the root, the
     * first 9 vertices of layer 2, the first 98 vertices of layer 3 or the
     * first 980 vertices of layer 4; the betweenness of these vertices should
     * all be zeros.
     *
     * Also, the betweenness of the common grand-grandparent in layer 2 is easy
     * to calculate as any shortest path going between grand-grandchildren
     * reachable via its left child and via its right child should pass through
     * it. This gives us a betweenness of 100 * 100 = 10000 for this node.
     * Similar calculations reveal that the betweenness of its two children
     * will be 100 * 100 + 100 * 90 / 2 = 14500 (the same paths as above plus any
     * path that goes from any of its subtree to any other subtree). These are
     * also the maximal betweennesses. The betwennesses in the remaining layers
     * are 1945 and 199 if the vertex being considered is a descendant of the
     * common grand-grandparent in layer 2, and zero otherwise. */

    igraph_vs_range(&vs_source, 10911, 11111);
    igraph_vs_range(&vs_target, 10911, 11111);
    igraph_vector_init(&bet, 0);

    igraph_betweenness_subset(
        /* graph=     */ &g,
        /* res=       */ &bet,
        /* vids=      */ igraph_vss_all(),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources =  */ vs_source,
        /* target =   */ vs_target,
        /* weights=   */ NULL
    );

    printf("Max betweenness: %f\n", igraph_vector_max(&bet));

    n = igraph_vcount(&g);
    for (i = 0; i < n; i++) {
        igraph_integer_t expected;

        if (i >= 10911) {
            /* layer 5, in the subset. There are 199 shortest paths that
             * contain these nodes, but the nodes are the _endpoints_ of the
             * paths so they don't count in the betweenness score */
            expected = 0;
        } else if (i >= 1111) {
            /* layer 5, not in the subset */
            expected = 0;
        } else if (i >= 1091) {
            /* layer 4, rightmost 20 nodes */
            expected = 1945;
        } else if (i >= 111) {
            /* layer 4, remaining nodes */
            expected = 0;
        } else if (i >= 109) {
            /* layer 3, rightmost 2 nodes */
            expected = 14500;
        } else if (i >= 11) {
            /* layer 3, remaining nodes */
            expected = 0;
        } else if (i == 10) {
            /* layer 2, rightmost node */
            expected = 10000;
        } else {
            expected = 0;
        }

        if (VECTOR(bet)[i] != expected) {
            printf(
                "Invalid betweenness for vertex %" IGRAPH_PRId ", expected %" IGRAPH_PRId ", got %" IGRAPH_PRId "\n",
                i, expected, (igraph_integer_t) VECTOR(bet)[i]
            );
            break;
        }
    }

    igraph_vector_init(&bet2, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1.0);

    igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ IGRAPH_UNDIRECTED,
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

    printf("\nNon-trivial weighted graph using subset algorithm\n");
    printf("==========================================================\n");
    igraph_vector_int_view(&edges, nontriv, sizeof(nontriv) / sizeof(nontriv[0]));
    igraph_create(&g, &edges, 0, /* directed= */ IGRAPH_UNDIRECTED);
    igraph_vector_view(&weights, nontriv_weights,
                       sizeof(nontriv_weights) / sizeof(nontriv_weights[0]));

    igraph_vector_init(&bet, 0);
    igraph_betweenness_subset(/* graph=     */ &g,
        /* res=       */ &bet,
        /* vids=      */ igraph_vss_all(),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ igraph_vss_all(),
        /* target = */ igraph_vss_all(),
        /* weights=   */ &weights);

    print_vector(&bet);

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

    for (i = 0; i < 5; i++)
    {
        igraph_vector_int_init_range(&node_vec, 0, 5);
        igraph_vector_int_remove(&node_vec, i);
        igraph_vs_vector(&vs, &node_vec);
        igraph_vector_int_init_range(&source_vec, 0, 5);
        igraph_vector_int_remove(&source_vec, i);
        igraph_vs_vector(&vs_source, &source_vec);
        printf("subset without %" IGRAPH_PRId "\n", i);
        printf("Unweighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ vs,
            /* directed = */ IGRAPH_UNDIRECTED,
            /* sources = */ vs_source,
            /* target = */ igraph_vss_all(),
            /* weights=   */ NULL);
        print_vector(&bet);

        printf("Weighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ vs,
            /* directed = */ IGRAPH_UNDIRECTED,
            /* sources = */ vs_source,
            /* target = */ igraph_vss_all(),
            /* weights */ &weights);
        print_vector(&bet2);
        printf("\n");

        IGRAPH_ASSERT(igraph_vector_all_e(&bet, &bet2));
        igraph_vs_destroy(&vs);
        igraph_vector_int_destroy(&node_vec);
        igraph_vs_destroy(&vs_source);
        igraph_vector_int_destroy(&source_vec);
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

    for (i = 0; i < 3; i++)
    {
        igraph_vector_int_init_range(&target_vec, 0, 4);
        igraph_vector_int_remove(&target_vec, i);
        igraph_vs_vector(&vs_target, &target_vec);
        printf("subset without %" IGRAPH_PRId "\n", i);
        printf("Unweighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ igraph_vss_all(),
            /* target = */ vs_target,
            /* weights */ NULL);
        print_vector(&bet);

        printf("Weighted\n");
        igraph_betweenness_subset(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* sources = */ igraph_vss_all(),
            /* target = */ vs_target,
            /* weights */ &weights);
        print_vector(&bet2);
        printf("\n");

        igraph_vs_destroy(&vs_target);
        igraph_vector_int_destroy(&target_vec);
    }

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    printf("\nEmpty graph\n");
    printf("==========================================================\n");
    igraph_empty(&g, 2, IGRAPH_UNDIRECTED);
    igraph_vector_init(&bet, 0);
    igraph_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &bet,
        /* vids=      */ igraph_vss_all(),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ igraph_vss_all(),
        /* target = */ igraph_vss_all(),
        /* weights=   */ NULL);

    print_vector(&bet);

    igraph_vector_destroy(&bet);
    igraph_destroy(&g);

    printf("\n37x37 grid graph\n");
    printf("==========================================================\n");

    {
        igraph_vector_int_t dims;

        igraph_vector_int_init(&dims, 2);
        VECTOR(dims)[0] = 37;
        VECTOR(dims)[1] = 37;

        igraph_square_lattice(&g, &dims, 1, IGRAPH_UNDIRECTED, /* mutual */ 0, /* periodic */ 0);

        igraph_vector_init(&bet, 0);
        igraph_vector_int_init_range(&target_vec, 0, igraph_vcount(&g));
        igraph_vector_int_remove(&target_vec, 0);
        igraph_vs_vector(&vs_target, &target_vec);
        igraph_vector_int_init_range(&source_vec, 0, igraph_vcount(&g));
        igraph_vector_int_remove(&source_vec, 0);
        igraph_vs_vector(&vs_source, &source_vec);

        igraph_betweenness_subset (/* graph=     */ &g,
        /* res=       */ &bet,
        /* vids=      */ igraph_vss_all(),
        /* directed = */ IGRAPH_UNDIRECTED,
        /* sources = */ vs_source,
        /* target = */ vs_target,
        /* weights=   */ NULL);;
        printf("Max betweenness: %f\n", igraph_vector_max(&bet));

        igraph_vector_destroy(&bet);
        igraph_destroy(&g);
        igraph_vector_int_destroy(&dims);
        igraph_vs_destroy(&vs_target);
        igraph_vector_int_destroy(&target_vec);
        igraph_vs_destroy(&vs_source);
        igraph_vector_int_destroy(&source_vec);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
