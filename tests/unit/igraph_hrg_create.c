/*
   IGraph library.
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

int main(void) {
    igraph_hrg_t hrg;
    igraph_t graph, new_graph;
    igraph_vector_t prob;

    igraph_rng_seed(igraph_rng_default(), 42);

    printf("Two leaf nodes which are always connected:\n");
    igraph_vector_init_real(&prob, 1, 1.0);
    igraph_hrg_init(&hrg, 0);
    igraph_small(&graph, 3, IGRAPH_DIRECTED, 0,1, 0,2, -1);
    /* the graph and the prob form the hrg dendrogram with probabilities. So you're not just using some graph */
    /* igraph_hrg_fit is used if you just have some graph that you'd want to do HRG analysis on */
    igraph_hrg_create(&hrg, &graph, &prob);
    igraph_hrg_game(&new_graph, &hrg);
    print_graph_canon(&new_graph);
    igraph_destroy(&graph);
    igraph_destroy(&new_graph);
    igraph_vector_destroy(&prob);

    printf("Four leaf nodes, one node connected to all others:\n");
    igraph_vector_init_real(&prob, 3, 1.0, 0.0, 0.0);
    igraph_small(&graph, 7, IGRAPH_DIRECTED, 0,3, 0,1, 1,4, 1,2, 2,5, 2,6, -1);
    igraph_hrg_create(&hrg, &graph, &prob);
    igraph_hrg_game(&new_graph, &hrg);
    print_graph_canon(&new_graph);

    igraph_destroy(&graph);
    igraph_destroy(&new_graph);
    igraph_vector_destroy(&prob);

    VERIFY_FINALLY_STACK();

    printf("Check error handling for wrong number of probabilities.\n");
    igraph_vector_init_real(&prob, 3, 1.0, 0.0, 0.0);
    igraph_small(&graph, 3, IGRAPH_DIRECTED, 0,1, 0,2, -1);
    CHECK_ERROR(igraph_hrg_create(&hrg, &graph, &prob), IGRAPH_EINVAL);

    igraph_destroy(&graph);
    igraph_vector_destroy(&prob);

    printf("Check error handling for graph with loop.\n");
    igraph_vector_init_real(&prob, 1, 1.0);
    igraph_small(&graph, 3, IGRAPH_DIRECTED, 0,0, 0,2, -1);
    CHECK_ERROR(igraph_hrg_create(&hrg, &graph, &prob), IGRAPH_EINVAL);
    igraph_destroy(&graph);

    printf("Check error handling for graph with no edges.\n");
    igraph_small(&graph, 3, IGRAPH_DIRECTED, -1);
    CHECK_ERROR(igraph_hrg_create(&hrg, &graph, &prob), IGRAPH_EINVAL);
    igraph_destroy(&graph);
    igraph_vector_destroy(&prob);
    igraph_hrg_destroy(&hrg);

    VERIFY_FINALLY_STACK();

    return 0;
}
