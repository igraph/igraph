/* -*- mode: C++ -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

#include "test_utilities.inc"

int main() {
    igraph_t graph;
    igraph_t full, tree;
    igraph_hrg_t hrg;
    igraph_t dendrogram;
    // int i, j;
    // igraph_vector_t neis;

    igraph_rng_seed(igraph_rng_default(), 42);

    // We need attributes
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_full(&full, 10, /*directed=*/ 0, /*loops=*/ 0);
    igraph_tree(&tree, 15, /*children=*/ 2, /*type=*/ IGRAPH_TREE_UNDIRECTED);
    igraph_disjoint_union(&graph, &full, &tree);
    igraph_add_edge(&graph, 0, 10);

    igraph_destroy(&full);
    igraph_destroy(&tree);

    // Fit
    igraph_hrg_init(&hrg, igraph_vcount(&graph));
    igraph_hrg_fit(&graph, &hrg, /*start=*/ 0, /*steps=*/ 0);

    // Create a graph from it
    igraph_hrg_dendrogram(&dendrogram, &hrg);

    // Print the tree, with labels
    // igraph_vector_init(&neis, 0);
    // for (i=0; i<igraph_vcount(&graph)-1; i++) {
    //   printf("Vertex # %2i, ", (int) (i+igraph_vcount(&graph)));
    //   igraph_neighbors(&dendrogram, &neis, i+igraph_vcount(&graph), IGRAPH_OUT);
    //   printf("left: # %2i, right: # %2i, ", (int) VECTOR(neis)[0],
    //       (int) VECTOR(neis)[1]);
    //   printf("prob: %6.2g\n",
    //       VAN(&dendrogram, "probability", i+igraph_vcount(&graph)));
    // }
    // igraph_vector_destroy(&neis);

    igraph_destroy(&dendrogram);
    igraph_hrg_destroy(&hrg);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}

