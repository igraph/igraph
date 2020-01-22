/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

void prepare_weights_vector(igraph_vector_t* weights, const igraph_t* graph) {
    int i, n = igraph_ecount(graph);
    igraph_vector_resize(weights, n);
    for (i = 0; i < n; i++) {
        VECTOR(*weights)[i] = i % 5;
    }
}

int main() {
    igraph_t graph;

    igraph_vector_t v;
    igraph_real_t edges[] = { 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
                              0, 10, 0, 11, 0, 12, 0, 13, 0, 17, 0, 19, 0, 21, 0, 31,
                              1, 2, 1, 3, 1, 7, 1, 13, 1, 17, 1, 19, 1, 21, 1, 30,
                              2, 3, 2, 7, 2, 27, 2, 28, 2, 32, 2, 9, 2, 8, 2, 13,
                              3, 7, 3, 12, 3, 13, 4, 6, 4, 10, 5, 6, 5, 10, 5, 16,
                              6, 16, 8, 30, 8, 32, 8, 33, 9, 33, 13, 33, 14, 32, 14, 33,
                              15, 32, 15, 33, 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
                              22, 32, 22, 33, 23, 25, 23, 27, 23, 32, 23, 33, 23, 29,
                              24, 25, 24, 27, 24, 31, 25, 31, 26, 29, 26, 33, 27, 33,
                              28, 31, 28, 33, 29, 32, 29, 33, 30, 32, 30, 33, 31, 32, 31, 33,
                              32, 33
                            };

    igraph_vector_t membership;
    igraph_vector_t weights;
    igraph_real_t modularity;
    igraph_bool_t simple;
    int retval;

    igraph_vector_view(&v, edges, sizeof(edges) / sizeof(double));
    igraph_create(&graph, &v, 0, IGRAPH_UNDIRECTED);

    igraph_vector_init(&weights, 0);

    igraph_is_simple(&graph, &simple);
    if (!simple) {
        return 1;
    }

    igraph_vector_init(&membership, 0);

    igraph_set_error_handler(&igraph_error_handler_printignore);

    /* Zachary karate club, unweighted */
    retval = igraph_community_optimal_modularity(&graph, &modularity,
             &membership, 0);
    if (retval == IGRAPH_UNIMPLEMENTED) {
        return 77;
    }
    if (fabs(modularity - 0.4197896) > 0.0000001) {
        return 2;
    }
    /* Zachary karate club, weighted */
    prepare_weights_vector(&weights, &graph);
    igraph_community_optimal_modularity(&graph, &modularity,
                                        &membership, &weights);
    if (fabs(modularity - 0.5115767) > 0.0000001) {
        return 4;
    }
    igraph_destroy(&graph);

    /* simple graph with loop edges, unweighted */
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 2, 2, -1);
    igraph_community_optimal_modularity(&graph, &modularity,
                                        &membership, 0);
    if (fabs(modularity - 0.28125) > 0.00001) {
        return 3;
    }
    /* simple graph with loop edges, weighted */
    prepare_weights_vector(&weights, &graph);
    igraph_community_optimal_modularity(&graph, &modularity,
                                        &membership, &weights);
    if (fabs(modularity - 0.36686) > 0.00001) {
        return 5;
    }
    igraph_destroy(&graph);

    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&weights);

    return 0;
}

