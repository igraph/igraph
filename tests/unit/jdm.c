/*
   IGraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

void print_and_destroy(igraph_t* g, igraph_matrix_t* jdm, igraph_vector_t* weights) {
      print_matrix(jdm);

      igraph_destroy(g);
      igraph_matrix_destroy(jdm);
      if (weights) {
          igraph_vector_destroy(weights);
      }
}

int main (void) {
    igraph_t g;
    igraph_matrix_t jdm;
    igraph_vector_t weights;

    printf("Graph with no vertices\n");
    igraph_small(&g, 0, false, -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, -1, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Simple, undirected graph\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, -1, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Simple, directed graph\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Undirected, self-loops, no multiedges\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Directed, self-loops, no multiedges\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Undirected multigraph, no self-loops\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Directed multigraph, no self-loops\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Undirected multigraph with self-loops\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Directed multigraph with self-loops\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm, NULL);

    // Weight tests
    printf("Weighted, simple, undirected graph\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_vector_init_range(&weights, -1,6);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, -1, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, simple, directed graph\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_vector_init_range(&weights, -2,8);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, undirected, self-loops, no multiedges\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_vector_init_range(&weights, 1,9);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, directed, self-loops, no multiedges\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_vector_init_range(&weights, 1,12);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, undirected multigraph, no self-loops\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_vector_init_range(&weights, 1,9);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, directed multigraph, no self-loops\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);
    igraph_vector_init_range(&weights, 1,12);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, undirected multigraph with self-loops\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_vector_init_range(&weights, 1,10);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    printf("Weighted, directed multigraph with self-loops\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);
    igraph_vector_init_range(&weights, 1,13);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, &weights);
    print_and_destroy(&g, &jdm, &weights);

    // dout din tests

    // Directed
    printf("Directed: dout is small, cropped JDM (3x3)\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 3, 3, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Directed: din is small, cropped JDM (4, 2)\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 4, 2, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Directed: Automatic resize, dout < 0 (4x2)\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, -1, 2, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Directed: Valid dout and din (5x5)\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 5, 5, NULL);
    print_and_destroy(&g, &jdm, NULL);

    // Undirected
    printf("Undirected: dout is small, cropped JDM (3x4)\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 3, 4, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Undirected: din is small, cropped JDM (4x3)\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 4, 3, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Undirected: Automatic resize, dout or din < 0 (4x4)\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, -1, -1, NULL);
    print_and_destroy(&g, &jdm, NULL);

    printf("Undirected: Valid dout and din (5x5)\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 5, 5, NULL);
    print_and_destroy(&g, &jdm, NULL);

    // Clean up
    igraph_destroy(&g);
    igraph_matrix_destroy(&jdm);
    igraph_vector_destroy(&weights);
    VERIFY_FINALLY_STACK();

    return 0;
}
