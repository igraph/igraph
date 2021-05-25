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
#include <stdio.h>

#include "test_utilities.inc"

void print_and_destroy(igraph_t *graph, igraph_bool_t unconn) {
    igraph_real_t res, unconn_pairs;

    igraph_average_path_length(graph, &res, &unconn_pairs, IGRAPH_DIRECTED, unconn);
    printf("Average shortest path length: ");
    print_real(stdout, res, "%g");
    printf("\nNo. of unconnected pairs: %g\n", unconn_pairs);

    igraph_destroy(graph);
}

int main() {
    igraph_t graph;

    printf("Null graph:\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    print_and_destroy(&graph, 1);

    printf("\nSingleton graph:\n");
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    print_and_destroy(&graph, 1);

    printf("\n2-vertex empty graph, unconn=true:\n");
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    print_and_destroy(&graph, 1);

    printf("\n2-vertex empty graph, unconn=false:\n");
    igraph_empty(&graph, 2, IGRAPH_UNDIRECTED);
    print_and_destroy(&graph, 0);

    printf("\nSmallest bifurcating directed tree, unconn=true:\n");
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,1, 0,2, -1);
    print_and_destroy(&graph, 1);

    printf("\nSmallest bifurcating directed tree, unconn=false:\n");
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,1, 0,2, -1);
    print_and_destroy(&graph, 1);

    printf("\nUndirected 11-cycle:\n");
    igraph_ring(&graph, 11, IGRAPH_UNDIRECTED, 0, 1);
    print_and_destroy(&graph, 1);

    printf("\nDirected 6-cycle:\n");
    igraph_ring(&graph, 11, IGRAPH_DIRECTED, 0, 1);
    print_and_destroy(&graph, 1);

    /* Result verified by szhorvat on 2021-03-04. */
    printf("\nDe Bruijn graph, n=3 and m=5:\n");
    igraph_de_bruijn(&graph, 3, 5);
    print_and_destroy(&graph, 1);

    VERIFY_FINALLY_STACK();

    return 0;
}
