/*
   igraph library.
   Copyright (C) 2006-2025  The igraph development team <igraph@igraph.org>

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

int test_no_multiple(void) {
    igraph_t g;
    igraph_bool_t simple, has_multi;

    /* Ensure that the test is deterministic */
    igraph_rng_seed(igraph_rng_default(), 137);

    /* null graph */

    igraph_erdos_renyi_game_gnm(&g, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_destroy(&g);

    /* singleton */

    igraph_erdos_renyi_game_gnm(&g, 1, 0, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    IGRAPH_ASSERT(igraph_vcount(&g) == 1);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_destroy(&g);

    /* singleton with loop */

    igraph_erdos_renyi_game_gnm(&g, 1, 1, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

    IGRAPH_ASSERT(igraph_vcount(&g) == 1);
    IGRAPH_ASSERT(igraph_ecount(&g) == 1);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 1, 1, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);

    IGRAPH_ASSERT(igraph_vcount(&g) == 1);
    IGRAPH_ASSERT(igraph_ecount(&g) == 1);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_destroy(&g);


    /* directed with loops */
    igraph_erdos_renyi_game_gnm(&g, 10, 10 * 10 - 1, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 10 - 1);
    IGRAPH_ASSERT(igraph_is_directed(&g));

    igraph_has_multiple(&g, &has_multi); IGRAPH_ASSERT(! has_multi);

    igraph_simplify(&g, /*remove_multiple=*/ false, /*remove_loops=*/ true, /*edge_comb=*/ NULL);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 || igraph_ecount(&g) == 10 * 9 - 1);

    igraph_destroy(&g);

    /* directed without loops */
    igraph_erdos_renyi_game_gnm(&g, 10, 10 * 9 - 1, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 - 1);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    /* undirected with loops */
    igraph_erdos_renyi_game_gnm(&g, 10, 10 * 11 / 2 - 1, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 11 / 2 - 1);
    IGRAPH_ASSERT(! igraph_is_directed(&g));

    igraph_has_multiple(&g, &has_multi); IGRAPH_ASSERT(! has_multi);

    igraph_simplify(&g, /*remove_multiple=*/ false, /*remove_loops=*/ true, /*edge_comb=*/ NULL);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 / 2 || igraph_ecount(&g) == 10 * 9 / 2 - 1);

    igraph_destroy(&g);

    /* undirected without loops */
    igraph_erdos_renyi_game_gnm(&g, 10, 10 * 9 / 2 - 1, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10 * 9 / 2 - 1);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);


    /* Create a couple of large graphs too */
    igraph_erdos_renyi_game_gnm(&g, 100000, 200000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 100000, 200000, IGRAPH_DIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED); IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 100000, 200000, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_has_multiple(&g, &has_multi); IGRAPH_ASSERT(! has_multi);
    igraph_simplify(&g, false, true, /*edge_comb=*/ NULL);  /* only remove loops */
    igraph_destroy(&g);

    igraph_erdos_renyi_game_gnm(&g, 100000, 200000, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&g) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&g) == 200000);
    igraph_has_multiple(&g, &has_multi); IGRAPH_ASSERT(! has_multi);
    igraph_destroy(&g);

    return 0;
}

int test_multiple(void) {
    igraph_t graph1, graph2;
    igraph_bool_t same, has_loop;

    igraph_rng_seed(igraph_rng_default(), 137); /* Makes test deterministic */

    /* null graph */

    igraph_erdos_renyi_game_gnm(&graph1, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 0);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 0);
    print_graph(&graph1);
    igraph_destroy(&graph1);

    /* null graph with more than zero edges is invalid */

    CHECK_ERROR(
            igraph_erdos_renyi_game_gnm(&graph1, 0, 1, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED),
        IGRAPH_EINVAL
    );

    /* singleton with no loops and more than zero edges is invalid */

    CHECK_ERROR(
            igraph_erdos_renyi_game_gnm(&graph1, 1, 1, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED),
        IGRAPH_EINVAL
    );

    /* singleton with multiple loops */

    igraph_erdos_renyi_game_gnm(&graph1, 1, 2, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED); /* undirected with loops */
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 1);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 2);
    print_graph(&graph1);
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm(&graph1, 1, 3, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED); /* undirected with loops */
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 1);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 3);
    print_graph(&graph1);
    igraph_destroy(&graph1);

    /* directed with loops */

    igraph_erdos_renyi_game_gnm(&graph1, 10, 560, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 560);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_simplify(&graph1, true, false, NULL);
    IGRAPH_ASSERT(igraph_ecount(&graph1) <= 100);
    igraph_destroy(&graph1);

    /* directed without loops */

    igraph_erdos_renyi_game_gnm(&graph1, 10, 890, IGRAPH_DIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 890);
    IGRAPH_ASSERT(igraph_is_directed(&graph1));
    igraph_has_loop(&graph1, &has_loop); IGRAPH_ASSERT(! has_loop);
    igraph_simplify(&graph1, true, false, NULL);
    IGRAPH_ASSERT(igraph_ecount(&graph1) <= 90);
    igraph_destroy(&graph1);

    /* undirected with loops */

    igraph_erdos_renyi_game_gnm(&graph1, 10, 540, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 540);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_simplify(&graph1, true, false, NULL);
    IGRAPH_ASSERT(igraph_ecount(&graph1) <= 55);
    igraph_destroy(&graph1);

    /* undirected without loops */

    igraph_erdos_renyi_game_gnm(&graph1, 10, 440, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 10);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 440);
    IGRAPH_ASSERT(!igraph_is_directed(&graph1));
    igraph_has_loop(&graph1, &has_loop); IGRAPH_ASSERT(! has_loop);
    igraph_copy(&graph2, &graph1);
    igraph_simplify(&graph1, true, true, NULL);
    igraph_simplify(&graph2, true, false, NULL);
    igraph_isomorphic(&graph1, &graph2, &same);
    IGRAPH_ASSERT(same);
    igraph_destroy(&graph1);
    igraph_destroy(&graph2);

    /* large graphs */

    igraph_erdos_renyi_game_gnm(&graph1, 1000, 200000, IGRAPH_UNDIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 1000);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 200000);
    igraph_has_loop(&graph1, &has_loop); IGRAPH_ASSERT(! has_loop);
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm(&graph1, 100000, 200000, IGRAPH_DIRECTED, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 200000);
    igraph_has_loop(&graph1, &has_loop); IGRAPH_ASSERT(! has_loop);
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm(&graph1, 100000, 200000, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 200000);
    igraph_destroy(&graph1);

    igraph_erdos_renyi_game_gnm(&graph1, 100000, 200000, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
    IGRAPH_ASSERT(igraph_vcount(&graph1) == 100000);
    IGRAPH_ASSERT(igraph_ecount(&graph1) == 200000);
    igraph_destroy(&graph1);

    return 0;
}

int test_iea(void) {
    igraph_t g;
    igraph_bool_t has_loop;

    /* null graph */
    igraph_iea_game(&g, 0, 0, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    IGRAPH_ASSERT(!igraph_is_directed(&g));
    igraph_destroy(&g);

    /* null graph with more than zero edges is invalid */
    CHECK_ERROR(
        igraph_iea_game(&g, 0, 1, IGRAPH_UNDIRECTED, IGRAPH_LOOPS),
        IGRAPH_EINVAL
    );

    /* singleton with no loops and more than zero edges is invalid */
    CHECK_ERROR(
        igraph_iea_game(&g, 1, 1, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS),
        IGRAPH_EINVAL
    );

    /* singleton with loop */
    igraph_iea_game(&g, 1, 1, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 1);
    IGRAPH_ASSERT(igraph_ecount(&g) == 1);
    IGRAPH_ASSERT(!igraph_is_directed(&g));
    igraph_destroy(&g);

    /* singleton with three loops */
    igraph_iea_game(&g, 1, 3, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 1);
    IGRAPH_ASSERT(igraph_ecount(&g) == 3);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_destroy(&g);

    /* larger undirected */
    igraph_iea_game(&g, 5, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 5);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10);
    IGRAPH_ASSERT(!igraph_is_directed(&g));
    igraph_has_loop(&g, &has_loop); IGRAPH_ASSERT(! has_loop);
    igraph_destroy(&g);

    igraph_iea_game(&g, 5, 10, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 5);
    IGRAPH_ASSERT(igraph_ecount(&g) == 10);
    IGRAPH_ASSERT(!igraph_is_directed(&g));
    igraph_destroy(&g);

    /* larger directed */
    igraph_iea_game(&g, 4, 11, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 4);
    IGRAPH_ASSERT(igraph_ecount(&g) == 11);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_has_loop(&g, &has_loop); IGRAPH_ASSERT(! has_loop);
    igraph_destroy(&g);

    igraph_iea_game(&g, 4, 11, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    IGRAPH_ASSERT(igraph_vcount(&g) == 4);
    IGRAPH_ASSERT(igraph_ecount(&g) == 11);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    igraph_destroy(&g);

    return 0;
}

int main(void) {
    RUN_TEST(test_no_multiple);
    RUN_TEST(test_multiple);
    RUN_TEST(test_iea);
}
