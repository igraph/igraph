/*
   igraph library.
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

void print_and_destroy(igraph_t *g, igraph_neimode_t mode, igraph_vector_t *weights) {
    igraph_vector_t v;
    igraph_real_t value;

    igraph_vector_init(&v, 0);
    igraph_eigenvector_centrality(g, &v, &value,
                                  mode,
                                  weights,
                                  /* options */ NULL);

    printf("Eigenvalue: %g\n", value);
    printf("Eigenvector:\n");
    igraph_vector_print(&v);
    printf("\n");

    igraph_destroy(g);
    igraph_vector_destroy(&v);
}

int main(void) {

    igraph_t g;
    igraph_vector_t weights;

    printf("Undirected graph with no vertices:\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, -1);
    print_and_destroy(&g, IGRAPH_ALL, NULL);

    printf("Directed graph with no vertices:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
    print_and_destroy(&g, IGRAPH_OUT, NULL);

    printf("Undirected graph with no edges:\n");
    igraph_small(&g, 5, IGRAPH_UNDIRECTED, -1);
    print_and_destroy(&g, IGRAPH_ALL, NULL);

    printf("Directed graph with no edges:\n");
    igraph_small(&g, 5, IGRAPH_DIRECTED, -1);
    print_and_destroy(&g, IGRAPH_IN, NULL);

    printf("Undirected full graph:\n");
    igraph_full(&g, 5, IGRAPH_UNDIRECTED, /* loops */ false);
    print_and_destroy(&g, IGRAPH_ALL, NULL);

    printf("Directed full graph:\n");
    igraph_full(&g, 5, IGRAPH_DIRECTED, /* loops */ false);
    print_and_destroy(&g, IGRAPH_OUT, NULL);

    printf("Undirected full graph with weights:\n");
    igraph_vector_init(&weights, 10);
    igraph_vector_fill(&weights, 1);
    igraph_full(&g, 5, IGRAPH_UNDIRECTED, /* loops */ false);
    print_and_destroy(&g, IGRAPH_ALL, &weights);
    igraph_vector_destroy(&weights);

    printf("Directed full graph with weights:\n");
    igraph_full(&g, 5, IGRAPH_DIRECTED, /* loops */ false);
    igraph_vector_init(&weights, 20);
    igraph_vector_fill(&weights, 1);
    print_and_destroy(&g, IGRAPH_OUT, &weights);
    igraph_vector_destroy(&weights);

    printf("Undirected star graph with weights:\n");
    igraph_vector_init(&weights, 4);
    igraph_vector_fill(&weights, 1);
    igraph_star(&g, 5, IGRAPH_STAR_UNDIRECTED, 0);
    print_and_destroy(&g, IGRAPH_ALL, &weights);
    igraph_vector_destroy(&weights);

    printf("Directed out-star graph with weights:\n");
    igraph_star(&g, 5, IGRAPH_STAR_OUT, 0);
    igraph_vector_init(&weights, 4);
    igraph_vector_fill(&weights, 1);
    print_and_destroy(&g, IGRAPH_OUT, &weights);
    igraph_vector_destroy(&weights);

    printf("Directed in-star graph with weights:\n");
    igraph_star(&g, 5, IGRAPH_STAR_IN, 0);
    igraph_vector_init(&weights, 4);
    igraph_vector_fill(&weights, 2);
    print_and_destroy(&g, IGRAPH_OUT, &weights);
    igraph_vector_destroy(&weights);

    printf("Small directed graph, OUT:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 1,3,
                 -1);
    print_and_destroy(&g, IGRAPH_OUT, NULL);

    printf("Small directed graph, IN:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 1,3,
                 -1);
    print_and_destroy(&g, IGRAPH_IN, NULL);

    printf("Small directed graph, ALL:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 1,3,
                 -1);
    print_and_destroy(&g, IGRAPH_ALL, NULL);

    printf("Small directed graph with weights, OUT:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 1,3,
                 -1);
    igraph_vector_init_range(&weights, 1, igraph_ecount(&g)+1);
    print_and_destroy(&g, IGRAPH_OUT, &weights);
    igraph_vector_destroy(&weights);

    printf("Small directed graph with weights, IN:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 1,3,
                 -1);
    igraph_vector_init_range(&weights, 1, igraph_ecount(&g)+1);
    print_and_destroy(&g, IGRAPH_IN, &weights);
    igraph_vector_destroy(&weights);

    printf("Small directed graph with weights, ALL:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 1,3,
                 -1);
    igraph_vector_init_range(&weights, 1, igraph_ecount(&g)+1);
    print_and_destroy(&g, IGRAPH_ALL, &weights);
    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();

    printf("Check handling of wrong number of weights.\n");
    igraph_vector_init(&weights, 2);
    igraph_vector_fill(&weights, 1);
    igraph_full(&g, 5, IGRAPH_DIRECTED, /* loops */ false);
    CHECK_ERROR(igraph_eigenvector_centrality(&g, NULL, NULL, /* mode */ IGRAPH_OUT,
                                  &weights,
                                  /* options */ NULL), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_eigenvector_centrality(&g, NULL, NULL, /* mode */ IGRAPH_ALL,
                                  &weights,
                                  /* options */ NULL), IGRAPH_EINVAL);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
