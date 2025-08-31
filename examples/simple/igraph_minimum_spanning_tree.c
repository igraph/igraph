/*
   igraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t graph;
    igraph_vector_t eb;
    igraph_vector_int_t edges;

    /* Initialize the library. */
    igraph_setup();

    /* Create the vector where the tree edges will be stored. */
    igraph_vector_int_init(&edges, 0);

    /* Create the Frucht graph */
    igraph_famous(&graph, "Frucht");

    /* Compute the edge betweenness. */
    igraph_vector_init(&eb, igraph_ecount(&graph));
    igraph_edge_betweenness(&graph, /*weights=*/ NULL, &eb, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_UNDIRECTED, false);

    /* Use Prim's algorithm to compute the edges that belong to the minimum weight
     * spanning tree, using edge betweenness values as edge weights. */
    igraph_minimum_spanning_tree(&graph, &edges, &eb, IGRAPH_MST_PRIM);
    printf("Minimum spanning tree edges:\n");
    igraph_vector_int_print(&edges);

    /* A maximum spanning tree can be computed by first negating the weights. */
    igraph_vector_scale(&eb, -1);

    /* Compute and output the edges that belong to the maximum weight spanning tree,
     * letting igraph automatically select the most suitable algorithm. */
    igraph_minimum_spanning_tree(&graph, &edges, &eb, IGRAPH_MST_AUTOMATIC);
    printf("\nMaximum spanning tree edges:\n");
    igraph_vector_int_print(&edges);

    igraph_real_t total_tree_weight = 0;
    igraph_int_t n = igraph_vector_int_size(&edges);
    for (igraph_int_t i=0; i < n; i++) {
        total_tree_weight += -VECTOR(eb)[ VECTOR(edges)[i] ];
    }
    printf("\nTotal maximum spanning tree weight: %g\n", total_tree_weight);

    /* Clean up */
    igraph_destroy(&graph);
    igraph_vector_destroy(&eb);
    igraph_vector_int_destroy(&edges);

    return 0;
}
