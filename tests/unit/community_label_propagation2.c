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

/* Test case for bug #1852 */

int main(void) {
    igraph_t graph;
    igraph_vector_int_t membership, initial_labels;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Undirected graph with unlabelled components */

    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
                 1,2, -1);

    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_init(&initial_labels, igraph_vcount(&graph));
    VECTOR(initial_labels)[0] = 1;
    VECTOR(initial_labels)[1] = -1;
    VECTOR(initial_labels)[2] = -1;
    VECTOR(initial_labels)[3] = -1;

    igraph_community_label_propagation(&graph, &membership, IGRAPH_ALL, NULL, &initial_labels, NULL);
    print_vector_int(&membership);

    igraph_destroy(&graph);

    /* Directed graph with unlabelled nodes not reachable from any labelled ones. */

    igraph_small(&graph, 8, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 3, 1, 2, 4, 4, 5, 5, 2, 4, 6,
                 -1);

    igraph_vector_int_resize(&initial_labels, igraph_vcount(&graph));
    igraph_vector_int_null(&initial_labels);
    VECTOR(initial_labels)[0] = -1;
    VECTOR(initial_labels)[1] = -1;
    VECTOR(initial_labels)[2] = 1;
    VECTOR(initial_labels)[3] = -1;
    VECTOR(initial_labels)[4] = 2;
    VECTOR(initial_labels)[6] = -1;
    VECTOR(initial_labels)[7] = -1;

    igraph_community_label_propagation(&graph, &membership, IGRAPH_OUT, NULL, &initial_labels, NULL);
    print_vector_int(&membership);

    igraph_destroy(&graph);

    /* None of the nodes are labelled initially */

    igraph_full(&graph, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_resize(&initial_labels, igraph_vcount(&graph));
    igraph_vector_int_fill(&initial_labels, -1);

    igraph_community_label_propagation(&graph, &membership, IGRAPH_OUT, NULL, &initial_labels, NULL);
    print_vector_int(&membership);

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&initial_labels);
    igraph_vector_int_destroy(&membership);

    VERIFY_FINALLY_STACK();

    return 0;
}
