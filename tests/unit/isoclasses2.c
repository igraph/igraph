/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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

/* Check that isoclass() and isoclass_create() are consistent with each other. */
void verify_classes(void) {
    igraph_integer_t class;
    igraph_integer_t size;

    igraph_integer_t classcountD[] = { 1, 1, 3, 16, 218 }; /* no. of unlabelled directed graphs */
    igraph_integer_t classcountU[] = { 1, 1, 2, 4, 11, 34, 156 }; /* no. of unlabelled undirected graphs */

    /* Directed */
    for (size=3; size <= 4; size++) {
        for (class=0; class < classcountD[size]; class++) {
            igraph_t g;
            igraph_integer_t class2;

            igraph_isoclass_create(&g, size, class, IGRAPH_DIRECTED);
            igraph_isoclass(&g, &class2);
            igraph_destroy(&g);

            IGRAPH_ASSERT(class == class2);
        }
    }

    /* Undirected */
    for (size=3; size <= 6; size++) {
        for (class=0; class < classcountU[size]; class++) {
            igraph_t g;
            igraph_integer_t class2;

            igraph_isoclass_create(&g, size, class, IGRAPH_UNDIRECTED);
            igraph_isoclass(&g, &class2);
            igraph_destroy(&g);

            IGRAPH_ASSERT(class == class2);
        }
    }
}

/* Generate small random graphs and check that their isoclasses are identified correctly. */
void random_test(void) {
    igraph_integer_t size, i;

    igraph_rng_seed(igraph_rng_default(), 137);

    /* Directed */
    for (size=3; size <= 4; size++) {
        for (i=0; i < 200; ++i) {
            igraph_t g1, g2;
            igraph_integer_t class;
            igraph_bool_t iso;

            igraph_erdos_renyi_game_gnp(&g1, size, 0.5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

            igraph_isoclass(&g1, &class);
            igraph_isoclass_create(&g2, size, class, IGRAPH_DIRECTED);

            igraph_isomorphic_bliss(&g1, &g2, NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
            IGRAPH_ASSERT(iso);

            igraph_destroy(&g2);
            igraph_destroy(&g1);
        }
    }

    /* Undirected */
    for (size=3; size <= 6; size++) {
        for (i=0; i < 200; ++i) {
            igraph_t g1, g2;
            igraph_integer_t class;
            igraph_bool_t iso;

            igraph_erdos_renyi_game_gnp(&g1, size, 0.5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

            igraph_isoclass(&g1, &class);
            igraph_isoclass_create(&g2, size, class, IGRAPH_UNDIRECTED);

            igraph_isomorphic_bliss(&g1, &g2, NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
            IGRAPH_ASSERT(iso);

            igraph_destroy(&g2);
            igraph_destroy(&g1);
        }
    }
}

/* Generate a random graph, select random subgraphs, and check that their
 * isoclasses are identified correctly. */
void random_subgraph_test(void) {
    igraph_t graph;
    igraph_integer_t size, i;
    igraph_vector_int_t vids;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init(&vids, 0);

    /* Directed */

    igraph_erdos_renyi_game_gnp(&graph, 40, 0.5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

    for (size=3; size <= 4; size++) {
        for (i=0; i < 100; ++i) {
            igraph_t sg1, sg2;
            igraph_integer_t class;
            igraph_bool_t iso;

            igraph_random_sample(&vids, 0, igraph_vcount(&graph) - 1, size);
            igraph_isoclass_subgraph(&graph, &vids, &class);
            igraph_isoclass_create(&sg1, size, class, igraph_is_directed(&graph));
            igraph_induced_subgraph(&graph, &sg2, igraph_vss_vector(&vids), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);

            igraph_isomorphic_bliss(&sg1, &sg2, NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
            IGRAPH_ASSERT(iso);

            igraph_destroy(&sg1);
            igraph_destroy(&sg2);
        }
    }

    igraph_destroy(&graph);

    /* Undirected */

    igraph_erdos_renyi_game_gnp(&graph, 60, 0.5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    for (size=3; size <= 6; size++) {
        for (i=0; i < 100; ++i) {
            igraph_t sg1, sg2;
            igraph_integer_t class;
            igraph_bool_t iso;

            igraph_random_sample(&vids, 0, igraph_vcount(&graph) - 1, size);
            igraph_isoclass_subgraph(&graph, &vids, &class);
            igraph_isoclass_create(&sg1, size, class, igraph_is_directed(&graph));
            igraph_induced_subgraph(&graph, &sg2, igraph_vss_vector(&vids), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);

            igraph_isomorphic_bliss(&sg1, &sg2, NULL, NULL, &iso, NULL, NULL, IGRAPH_BLISS_F, NULL, NULL);
            IGRAPH_ASSERT(iso);

            igraph_destroy(&sg1);
            igraph_destroy(&sg2);
        }
    }

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&vids);
}

int main(void) {

    verify_classes();
    random_test();
    random_subgraph_test();

    VERIFY_FINALLY_STACK();

    return 0;
}
