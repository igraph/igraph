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

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("No vertices:\n");
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/0, /*power*/ 0.0,
                  /*window*/ 1, /*edges per step(m)*/ 1, /*outseq*/ NULL, /*outpref?*/ 0, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("No edges:\n");
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/ 5, /*power*/ 0.0,
                  /*window*/ 1, /*edges per step(m)*/ 0, /*outseq*/ NULL, /*outpref?*/ 1, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    printf("A star with double edges.\n");
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/ 10, /*power*/ 30.0,
                  /*window*/ 100, /*edges per step(m)*/ 2, /*outseq*/ NULL, /*outpref?*/ 0, /*zero appeal*/ 0.001,
                  /*directed?*/ 1) == IGRAPH_SUCCESS);
    print_graph_canon(&g);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    /*Negative number of vertices*/
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/ -1, /*power*/ 10.0,
                  /*window*/ 100, /*edges per step(m)*/ 1, /*outseq*/ NULL, /*outpref?*/ 0, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    /*Negative number of edges*/
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/ 1, /*power*/ 10.0,
                  /*window*/ 100, /*edges per step(m)*/ -1, /*outseq*/ NULL, /*outpref?*/ 0, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    /*Negative window*/
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/ 1, /*power*/ 10.0,
                  /*window*/ -100, /*edges per step(m)*/ 1, /*outseq*/ NULL, /*outpref?*/ 0, /*zero appeal*/ 1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    /*Negative zero appeal*/
    IGRAPH_ASSERT(igraph_recent_degree_game(&g, /*number of vertices (n)*/ 1, /*power*/ 10.0,
                  /*window*/ 100, /*edges per step(m)*/ 1, /*outseq*/ NULL, /*outpref?*/ 0, /*zero appeal*/ -1,
                  /*directed?*/ 0) == IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();
    return 0;
}
