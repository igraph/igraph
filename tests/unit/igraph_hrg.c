/* -*- mode: C++ -*-  */
/*
   igraph library.
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

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_t full, tree;
    igraph_hrg_t hrg;
    igraph_t dendrogram;
    igraph_vector_t prob;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_full(&full, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_kary_tree(&tree, 15, /*children=*/ 2, /*type=*/ IGRAPH_TREE_UNDIRECTED);
    igraph_disjoint_union(&graph, &full, &tree);
    igraph_add_edge(&graph, 0, 10);

    igraph_destroy(&full);
    igraph_destroy(&tree);

    // Fit
    igraph_hrg_init(&hrg, igraph_vcount(&graph));
    igraph_hrg_fit(&graph, &hrg, /*start=*/ false, /*steps=*/ 0);

    // Create a graph from it
    igraph_vector_init(&prob, 0);
    igraph_from_hrg_dendrogram(&dendrogram, &hrg, &prob);

    // Print the tree, with labels
    igraph_vector_int_t neis;
    igraph_vector_int_init(&neis, 0);
    for (igraph_int_t i=0; i < igraph_vcount(&graph)-1; i++) {
        printf("Vertex # %2" IGRAPH_PRId ", ", (i+igraph_vcount(&graph)));
        igraph_neighbors(
            &dendrogram, &neis, i+igraph_vcount(&graph), IGRAPH_OUT,
            IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE
        );
        printf("left: # %2" IGRAPH_PRId ", right: # %2" IGRAPH_PRId ", ", VECTOR(neis)[0], VECTOR(neis)[1]);
        printf("prob: %6.2g\n", VECTOR(prob)[i+igraph_vcount(&graph)]);
    }
    igraph_vector_int_destroy(&neis);

    igraph_vector_destroy(&prob);
    igraph_destroy(&dendrogram);
    igraph_hrg_destroy(&hrg);
    igraph_destroy(&graph);

    // test small graph, omit probabilities
    igraph_small(&graph, 3, IGRAPH_UNDIRECTED,
                 0,1,
                 -1);

    igraph_hrg_init(&hrg, igraph_vcount(&graph));
    igraph_hrg_fit(&graph, &hrg, /*start=*/ false, /*steps=*/ 0);
    igraph_from_hrg_dendrogram(&dendrogram, &hrg, NULL);
    igraph_destroy(&dendrogram);
    igraph_hrg_destroy(&hrg);

    // test with a specific number of steps
    igraph_hrg_init(&hrg, igraph_vcount(&graph));
    igraph_hrg_fit(&graph, &hrg, /*start=*/ false, /*steps=*/ 10);
    igraph_from_hrg_dendrogram(&dendrogram, &hrg, NULL);
    igraph_destroy(&dendrogram);
    igraph_hrg_destroy(&hrg);

    igraph_destroy(&graph);

    // graph must have at least 3 vertices at the moment
    igraph_full(&graph, 2, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_hrg_init(&hrg, igraph_vcount(&graph));
    CHECK_ERROR(igraph_hrg_fit(&graph, &hrg, /*start=*/ false, /*steps=*/ 0), IGRAPH_EINVAL);
    igraph_hrg_destroy(&hrg);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
