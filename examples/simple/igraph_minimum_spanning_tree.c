/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int main(void) {

    igraph_t graph, tree;
    igraph_vector_t eb;
    igraph_vector_int_t edges;

    /* Create the Frucht graph */
    igraph_famous(&graph, "Frucht");

    /* Compute the edge betweenness. */
    igraph_vector_init(&eb, igraph_ecount(&graph));
    igraph_edge_betweenness(&graph, &eb, IGRAPH_UNDIRECTED, /*weights=*/ NULL);

    /* Compute and output a minimum weight spanning tree using edge betweenness
     * values as weights. */
    igraph_minimum_spanning_tree_prim(&graph, &tree, &eb);
    printf("Minimum spanning tree:\n");
    igraph_write_graph_edgelist(&tree, stdout);

    /* A maximum spanning tree can be computed by first negating the weights. */
    igraph_vector_scale(&eb, -1);

    /* Compute and output the edges that belong to the maximum weight spanning tree. */
    igraph_vector_int_init(&edges, 0);
    igraph_minimum_spanning_tree(&graph, &edges, &eb);
    printf("\nMaximum spanning tree edges:\n");
    igraph_vector_int_print(&edges);

    igraph_real_t total_tree_weight = 0;
    igraph_integer_t n = igraph_vector_int_size(&edges);
    for (igraph_integer_t i=0; i < n; i++) {
        total_tree_weight += -VECTOR(eb)[ VECTOR(edges)[i] ];
    }
    printf("\nTotal maximum spanning tree weight: %g\n", total_tree_weight);

    /* Clean up */
    igraph_vector_int_destroy(&edges);
    igraph_destroy(&tree);
    igraph_destroy(&graph);
    igraph_vector_destroy(&eb);

    return 0;
}
