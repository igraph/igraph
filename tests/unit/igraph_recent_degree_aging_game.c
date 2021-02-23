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
#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_bool_t tree;
    igraph_vector_t outseq;
    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/0, /*edges per step(m)*/ 1, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 1, /*aging_exp*/ 1, /*aging_bins*/ 1, /*time_window*/ 1, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("No edges:\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/5, /*edges per step(m)*/ 0, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 1, /*aging_exp*/ 1, /*aging_bins*/ 6, /*time_window*/ 1, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Prefer more edges to make a star of double edges:\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/5, /*edges per step(m)*/ 2, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 20, /*aging_exp*/ 0, /*aging_bins*/ 1, /*time_window*/ 100, /*zero appeal*/ 0.001,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Prefer older edges to make a star:\n");
    igraph_vector_init_int(&outseq, 7, 1, 2, 1, 2, 1, 2, 1);
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/7, /*edges per step(m)*/ 0, &outseq,
                  /*outpref?*/ 0, /*pa_exp*/ 0, /*aging_exp*/ 20, /*aging_bins*/ 8, /*time_window*/ 100, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("Checking if adding one edge per step makes a tree.\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/10, /*edges per step(m)*/ 1, /*outseq*/ NULL,
                  /*outpref?*/ 1, /*pa_exp*/ 2, /*aging_exp*/ 2, /*aging_bins*/ 5, /*time_window*/ 4, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    igraph_is_tree(&g, &tree, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(tree);
    igraph_destroy(&g);
    igraph_vector_destroy(&outseq);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking if these errors are properly handled:\n");
    printf("Negative number of vertices.\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/-20, /*edges per step(m)*/ 1, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 2, /*aging_exp*/ 2, /*aging_bins*/ 3, /*time_window*/ 4, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    printf("Negative number of edges.\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/20, /*edges per step(m)*/ -1, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 2, /*aging_exp*/ 2, /*aging_bins*/ 10, /*time_window*/ 4, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    printf("Negative aging bin.\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/20, /*edges per step(m)*/ 1, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 2, /*aging_exp*/ 2, /*aging_bins*/ -3, /*time_window*/ 4, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    printf("Negative zero appeal.\n");
    IGRAPH_ASSERT(igraph_recent_degree_aging_game(&g, /*nodes*/20, /*edges per step(m)*/ 1, /*outseq*/ NULL,
                  /*outpref?*/ 0, /*pa_exp*/ 2, /*aging_exp*/ 2, /*aging_bins*/ 10, /*time_window*/ 4, /*zero appeal*/ -1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);


    VERIFY_FINALLY_STACK();
    return 0;
}
