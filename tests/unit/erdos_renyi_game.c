/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_bool_t simple;

    /* Ensure that the test is deterministic */
    igraph_rng_seed(igraph_rng_default(), 137);

    /* G(n,p) */

    /* Empty graph */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.0,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.0,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    /* Complete graph */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 / 2);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 11 / 2);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(! simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 1.0,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 10);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(! simple);
    igraph_destroy(&g);

    /* Random graph */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.5,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.5,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.5,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 10, 0.5,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_destroy(&g);

    /* Create a couple of large graphs too */

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100000, 2.0 / 100000,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    igraph_destroy(&g);


    /* --------------------------------------------------------------------- */
    /* G(n,m) */

    /* directed with loops */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 10 - 1,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 10 - 1);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ NULL);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 || igraph_ecount(&g) == 10 * 9 - 1);

    igraph_destroy(&g);

    /* directed without loops */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 9 - 1,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 - 1);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    /* undirected with loops */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 11 / 2 - 1,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 11 / 2 - 1);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_simplify(&g, /*multiple=*/0, /*loops=*/1, /*edge_comb=*/ NULL);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 / 2 || igraph_ecount(&g) == 10 * 9 / 2 - 1);

    igraph_destroy(&g);

    /* undirected without loops */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 10, 10 * 9 / 2 - 1,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 / 2 - 1);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);


    /* Create a couple of large graphs too */
    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_is_simple(&g, &simple); IGRAPH_ASSERT(simple);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_simplify(&g, 0, 1, /*edge_comb=*/ 0);  /* only remove loops */
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 100000, 2.0 * 100000,
                            IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
