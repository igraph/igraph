/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

void print_check_destroy(igraph_t *graph, igraph_vector_int_list_t *result) {
    igraph_int_t i, rank = igraph_vector_int_list_size(result);
    igraph_int_t ecount = igraph_ecount(graph), vcount = igraph_vcount(graph);
    igraph_int_t ccount;

    for (i=0; i < rank; ++i) {
        igraph_vector_int_t *cycle = igraph_vector_int_list_get_ptr(result, i);
        IGRAPH_ASSERT(igraph_vector_int_size(cycle) > 0);
        igraph_int_t first_edge = VECTOR(*cycle)[0];
        igraph_int_t last_edge = VECTOR(*cycle)[igraph_vector_int_size(cycle) - 1];
        /* Verify that first and last edges of the cycle share a vertex */
        IGRAPH_ASSERT(
            (IGRAPH_FROM(graph, first_edge) == IGRAPH_FROM(graph, last_edge) ||
             IGRAPH_FROM(graph, first_edge) == IGRAPH_TO(graph, last_edge)) ||
            (IGRAPH_TO(graph, first_edge) == IGRAPH_FROM(graph, last_edge) ||
             IGRAPH_TO(graph, first_edge) == IGRAPH_TO(graph, last_edge))
        );
        print_vector_int(cycle);
    }
    igraph_connected_components(graph, NULL, NULL, &ccount, IGRAPH_WEAK);
    IGRAPH_ASSERT(rank == ecount - vcount + ccount);
    igraph_destroy(graph);
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t result;
    igraph_int_t rank;

    igraph_vector_int_list_init(&result, 0);

    printf("FUNDAMENTAL CYCLE BASIS\n");

    printf("\nNull graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_fundamental_cycles(&graph, NULL, &result, -1, -1);
    print_check_destroy(&graph, &result);

    printf("\nSingleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_fundamental_cycles(&graph, NULL, &result, -1, -1);
    print_check_destroy(&graph, &result);

    printf("\nSingle vertex with loop\n");
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0,0,
                 -1);
    igraph_fundamental_cycles(&graph, NULL, &result, -1, -1);
    print_check_destroy(&graph, &result);

    printf("\nTree\n");
    igraph_kary_tree(&graph, 3, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_fundamental_cycles(&graph, NULL, &result, -1, -1);
    print_check_destroy(&graph, &result);

    printf("\n2-cycle\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 0,1,
                 -1);
    igraph_fundamental_cycles(&graph, NULL, &result, -1, -1);
    print_check_destroy(&graph, &result);

    printf("\nDisconnected\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 1,2, 2,3, 3,1,
                 4,5, 5,4, 4,5,
                 6,7, 7,8, 8,9, 9,6, 6,8,
                 10,10, 10,11,
                 12,12,
                 -1);
    igraph_fundamental_cycles(&graph, NULL, &result, -1, -1);
    print_check_destroy(&graph, &result);

    printf("\nMINIMUM WEIGHT CYCLE BASIS\n");

    printf("\nNull graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);
    print_check_destroy(&graph, &result);

    printf("\nSingleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);
    print_check_destroy(&graph, &result);

    printf("\nSingle vertex with loop\n");
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED,
                 0,0,
                 -1);
    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);
    print_check_destroy(&graph, &result);

    printf("\nTree\n");
    igraph_kary_tree(&graph, 3, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);
    print_check_destroy(&graph, &result);

    printf("\n2-cycle\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 0,1,
                 -1);
    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);
    print_check_destroy(&graph, &result);

    printf("\nDisconnected\n");
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 1,2, 2,3, 3,1,
                 4,5, 5,4, 4,5,
                 6,7, 7,8, 8,9, 9,6, 6,8,
                 10,10, 10,11,
                 12,12,
                 -1);
    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);
    print_check_destroy(&graph, &result);

    printf("\nPeriodic (5,6)-grid\n");


    {
        igraph_vector_int_t dimvec;
        igraph_vector_bool_t periodic;

        igraph_vector_int_init_int_end(&dimvec, -1,
                                       1, 5, 6, -1);

        igraph_vector_bool_init(&periodic, igraph_vector_int_size(&dimvec));
        igraph_vector_bool_fill(&periodic, 1);

        igraph_square_lattice(&graph, &dimvec, 1, IGRAPH_ADJ_UNDIRECTED, /* mutual */ false, &periodic);

        igraph_vector_bool_destroy(&periodic);
        igraph_vector_int_destroy(&dimvec);
    }

    igraph_minimum_cycle_basis(&graph, NULL, &result, /* cutoff */ -1, /* complete */ true, /* ordered */ true);

    rank = igraph_vector_int_list_size(&result);

    /* In a periodic grid graph, all elements in the minimum basis have size 4,
     * except two, which have size equal to the grid dimensions. */

    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[0]) == 4);
    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[rank-1]) == 6);
    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[rank-2]) == 5);
    IGRAPH_ASSERT(igraph_vector_int_size(&VECTOR(result)[rank-3]) == 4);

    print_check_destroy(&graph, &result);

#if 0
    /* This is for benchmarking */
    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_watts_strogatz_game(&graph, 3, 10, 2, 0.1, 0, 0);
    igraph_minimum_cycle_basis(&graph, /* cutoff */ -1, /* complete */ true, * ordered */ true, &result);
    print_check_destroy(&graph, &result);
#endif

    igraph_vector_int_list_destroy(&result);

    VERIFY_FINALLY_STACK();

    return 0;
}
