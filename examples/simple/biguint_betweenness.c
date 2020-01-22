/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

int check(const igraph_vector_t *v1, const igraph_vector_t *v2, int code) {
    igraph_vector_t v;
    long int i, n = igraph_vector_size(v1);
    igraph_real_t m;

    igraph_vector_copy(&v, v1);
    igraph_vector_sub(&v, v2);

    for (i = 0; i < n; i++) {
        VECTOR(v)[i] = fabs(VECTOR(v)[i]);
    }

    if ( (m = igraph_vector_max(&v)) > 0.01) {
        printf("Difference: %g\n", m);
        exit(code);
    }

    igraph_vector_destroy(&v);

    return 0;
}

int main() {

    igraph_t g;
    igraph_vector_t bet, bet2, weights, edges;
    igraph_vector_t bbet, bbet2;

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
    igraph_vector_init(&bbet, 0);

    igraph_betweenness_estimate(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* cutoff=    */ 2,
            /* weights=   */ 0,
            /* nobigint=  */ 1);

    igraph_betweenness_estimate(/* graph=     */ &g,
            /* res=       */ &bbet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* cutoff=    */ 2,
            /* weights=   */ 0,
            /* nobigint=  */ 0);

    check(&bet, &bbet, 10);

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bbet);
    igraph_destroy(&g);

    /*******************************************************/

    igraph_tree(&g, 20000, 10, IGRAPH_TREE_UNDIRECTED);

    igraph_vector_init(&bet, 0);
    igraph_vector_init(&bbet, 0);

    igraph_betweenness_estimate(/* graph=     */ &g,
            /* res=       */ &bet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* cutoff=    */ 3,
            /* weights=   */ 0,
            /* nobigint=  */ 1);

    igraph_betweenness_estimate(/* graph=     */ &g,
            /* res=       */ &bbet,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* cutoff=    */ 3,
            /* weights=   */ 0,
            /* nobigint=  */ 0);

    check(&bet, &bbet, 20);

    igraph_vector_init(&bet2, 0);
    igraph_vector_init(&bbet2, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1.0);

    igraph_betweenness_estimate(/* graph=     */ &g,
            /* res=       */ &bet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* cutoff=    */ 3,
            /* weights=   */ &weights,
            /* nobigint=  */ 1);

    igraph_betweenness_estimate(/* graph=     */ &g,
            /* res=       */ &bbet2,
            /* vids=      */ igraph_vss_all(),
            /* directed = */ 0,
            /* cutoff=    */ 3,
            /* weights=   */ &weights,
            /* nobigint=  */ 0);

    if (!igraph_vector_all_e(&bet, &bet2)) {
        return 1;
    }

    /*   if (!igraph_vector_all_e(&bbet, &bbet2)) { */
    /*     return 2; */
    /*   } */

    check(&bet, &bbet, 30);
    check(&bet2, &bbet2, 40);

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bet2);
    igraph_vector_destroy(&bbet);
    igraph_vector_destroy(&bbet2);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    /* Non-trivial weighted graph */
    igraph_vector_view(&edges, nontriv, sizeof(nontriv) / sizeof(igraph_real_t));
    igraph_create(&g, &edges, 0, /* directed= */ 0);
    igraph_vector_view(&weights, nontriv_weights,
                       sizeof(nontriv_weights) / sizeof(igraph_real_t));
    igraph_vector_init(&bet, 0);
    igraph_vector_init(&bbet, 0);

    igraph_betweenness(/*graph=*/ &g, /*res=*/ &bet, /*vids=*/ igraph_vss_all(),
                                  /*directed=*/0, /*weights=*/ &weights, /*nobigint=*/ 1);

    igraph_betweenness(/*graph=*/ &g, /*res=*/ &bbet, /*vids=*/ igraph_vss_all(),
                                  /*directed=*/0, /*weights=*/ &weights, /*nobigint=*/ 0);

    igraph_vector_view(&bet2, nontriv_res,
                       sizeof(nontriv_res) / sizeof(igraph_real_t));

    if (!igraph_vector_all_e(&bet, &bet2)) {
        return 2;
    }

    check(&bet, &bbet, 50);

    igraph_vector_destroy(&bet);
    igraph_vector_destroy(&bbet);
    igraph_destroy(&g);

    if (IGRAPH_FINALLY_STACK_SIZE() != 0) {
        return 3;
    }

    return 0;
}
