/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

int main() {

    igraph_t g;
    int i;
    igraph_bool_t simple;

    /* G(n,p) */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.0,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    if (igraph_ecount(&g) != 0) {
        return 1;
    }
    if (igraph_is_directed(&g)) {
        return 2;
    }
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    if (igraph_ecount(&g) != 10 * 9) {
        return 3;
    }
    if (!igraph_is_directed(&g)) {
        return 4;
    }
    igraph_destroy(&g);

    /* More useful tests */
    /*   printf("directed with loops\n"); */
    for (i = 0; i < 100; i++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.9999999,
                                IGRAPH_DIRECTED, IGRAPH_LOOPS);
        if (igraph_vcount(&g) != 10) {
            return 5;
        }
        if (igraph_ecount(&g) != 10 * 10) {
            return 77;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        if (igraph_ecount(&g) != 10 * 9) {
            return 77;
        }
        igraph_destroy(&g);
    }

    /*   printf("directed without loops\n"); */
    for (i = 0; i < 100; i++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.9999999,
                                IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
        if (igraph_vcount(&g) != 10) {
            return 7;
        }
        if (igraph_ecount(&g) != 10 * (10 - 1)) {
            return 77;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        if (igraph_ecount(&g) != 10 * 9) {
            return 77;
        }
        igraph_destroy(&g);
    }

    /*   printf("undirected with loops\n"); */
    for (i = 0; i < 100; i++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.9999999,
                                IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
        if (igraph_vcount(&g) != 10) {
            return 9;
        }
        if (igraph_ecount(&g) != 10 * (10 + 1) / 2) {
            return 77;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        if (igraph_ecount(&g) != 10 * (10 - 1) / 2) {
            return 77;
        }
        igraph_destroy(&g);
    }

    /*   printf("undirected without loops\n"); */
    for (i = 0; i < 100; i++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.9999999,
                                IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
        if (igraph_vcount(&g) != 10) {
            return 11;
        }
        if (igraph_ecount(&g) != 10 * (10 - 1) / 2) {
            return 77;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        if (igraph_ecount(&g) != 10 * (10 - 1) / 2) {
            return 77;
        }
        igraph_destroy(&g);
    }

    /* Create a couple of large graphs too */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 25;
    }
    igraph_destroy(&g);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 25;
    }
    igraph_destroy(&g);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 25;
    }
    igraph_destroy(&g);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 25;
    }
    igraph_destroy(&g);


    /* --------------------------------------------------------------------- */
    /* G(n,m) */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 0.5,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_destroy(&g);

    /* More useful tests */
    /*   printf("directed with loops\n"); */
    for (i = 0; i < 100; i++) {
        long int ec;
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 10 - 1,
                                IGRAPH_DIRECTED, IGRAPH_LOOPS);
        if (igraph_vcount(&g) != 10) {
            return 13;
        }
        if (igraph_ecount(&g) != 10 * 10 - 1) {
            return 14;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        igraph_is_simple(&g, &simple);
        if (!simple) {
            return 27;
        }
        ec = igraph_ecount(&g);
        if (ec != 10 * 9 && ec != 10 * 9 - 1) {
            return 15;
        }
        igraph_destroy(&g);
    }

    /*   printf("directed without loops\n"); */
    for (i = 0; i < 100; i++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 9 - 1,
                                IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
        igraph_is_simple(&g, &simple);
        if (!simple) {
            return 28;
        }
        if (igraph_vcount(&g) != 10) {
            return 16;
        }
        if (igraph_ecount(&g) != 10 * (10 - 1) - 1) {
            return 17;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        if (igraph_ecount(&g) != 10 * 9 - 1) {
            return 18;
        }
        igraph_destroy(&g);
    }

    /*   printf("undirected with loops\n"); */
    for (i = 0; i < 100; i++) {
        long int ec;
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 11 / 2 - 1,
                                IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
        if (igraph_vcount(&g) != 10) {
            return 19;
        }
        if (igraph_ecount(&g) != 10 * (10 + 1) / 2 - 1) {
            return 20;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        igraph_is_simple(&g, &simple);
        if (!simple) {
            return 29;
        }
        ec = igraph_ecount(&g);
        if (ec != 10 * (10 - 1) / 2 && ec != 10 * 9 / 2 - 1) {
            return 21;
        }
        igraph_destroy(&g);
    }

    /*   printf("undirected without loops\n"); */
    for (i = 0; i < 100; i++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 9 / 2 - 1,
                                IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
        igraph_is_simple(&g, &simple);
        if (!simple) {
            return 30;
        }
        if (igraph_vcount(&g) != 10) {
            return 22;
        }
        if (igraph_ecount(&g) != 10 * (10 - 1) / 2 - 1) {
            return 23;
        }
        igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ 0);
        if (igraph_ecount(&g) != 10 * (10 - 1) / 2 - 1) {
            return 24;
        }
        igraph_destroy(&g);
    }

    /* Create a couple of large graphs too */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 26;
    }
    if (igraph_ecount(&g) != 200000) {
        return 26;
    }
    igraph_is_simple(&g, &simple);
    if (!simple) {
        return 31;
    }
    igraph_destroy(&g);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_is_simple(&g, &simple);
    if (!simple) {
        return 32;
    }
    if (igraph_vcount(&g) != 100000) {
        return 26;
    }
    if (igraph_ecount(&g) != 200000) {
        return 26;
    }
    igraph_destroy(&g);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 26;
    }
    if (igraph_ecount(&g) != 200000) {
        return 26;
    }
    igraph_simplify(&g, 0, 1, /*edge_comb=*/ 0);  /* only remove loops */
    igraph_is_simple(&g, &simple);
    if (!simple) {
        return 33;
    }
    igraph_destroy(&g);
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    if (igraph_vcount(&g) != 100000) {
        return 26;
    }
    if (igraph_ecount(&g) != 200000) {
        return 26;
    }
    igraph_simplify(&g, 0, 1, /*edge_comb=*/ 0);  /* only remove loops */
    igraph_is_simple(&g, &simple);
    if (!simple) {
        return 34;
    }
    igraph_destroy(&g);

    return 0;
}
