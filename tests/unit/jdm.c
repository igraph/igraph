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

void print_and_destroy(igraph_t* g, igraph_matrix_int_t* jdm) {
      igraph_matrix_int_print(jdm);
      printf("\n");

      igraph_destroy(g);
      igraph_matrix_int_destroy(jdm);
}

int main (void) {
    // TODO: Structure tests and output file
    // Structures that need to be destroyed
    igraph_t g;
    igraph_matrix_int_t jdm;
    igraph_vector_int_t weights;

    printf("Graph with no vertices\n");
    igraph_small(&g, 0, false, -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, -1, 0, NULL);
//    igraph_matrix_int_print(&jdm);
    print_and_destroy(&g, &jdm);

    printf("Simple, undirected graph\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, -1, NULL);
    print_and_destroy(&g, &jdm);

    printf("Simple, directed graph\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    printf("Undirected, self-loops, no multiedges\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    printf("Directed, self-loops, no multiedges\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    printf("Undirected multigraph, no self-loops\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    printf("Directed multigraph, no self-loops\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    printf("Undirected multigraph with self-loops\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 4, 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    printf("Directed multigraph with self-loops\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 1, 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 0, 0, NULL);
    print_and_destroy(&g, &jdm);

    // TODO: Weights tests
    // Initialize weight vector
    igraph_vector_int_init(&weights, 0);
    igraph_vector_int_print(&weights);


    // dout din tests
    // Directed
    // Erroneous calls

    // dout is too small
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    CHECK_ERROR(igraph_construct_jdm(&g, &jdm, 3, 3, NULL), IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_matrix_int_destroy(&jdm);
    // din is too small
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    CHECK_ERROR(igraph_construct_jdm(&g, &jdm, 4, 2, NULL), IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_matrix_int_destroy(&jdm);

    // Successful calls
    printf("Automatic resize, dout or din < 0\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, -1, 2, NULL);
    print_and_destroy(&g, &jdm);

    printf("Valid dout and din\n");
    igraph_small(&g, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2, 4, 3,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 5, 5, NULL);
    print_and_destroy(&g, &jdm);

    // Undirected
    // Erroneous calls

    // dout is too small -> Throw error
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    CHECK_ERROR(igraph_construct_jdm(&g, &jdm, 3, 4, NULL), IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_matrix_int_destroy(&jdm);
    // din is too small -> Throw error
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    CHECK_ERROR(igraph_construct_jdm(&g, &jdm, 4, 3, NULL), IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_matrix_int_destroy(&jdm);

    // Successful calls
    printf("Automatic resize, dout or din < 0\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, -1, -1, NULL);
    print_and_destroy(&g, &jdm);

    printf("Valid dout and din\n");
    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);
    igraph_matrix_int_init(&jdm, 1, 1);
    igraph_construct_jdm(&g, &jdm, 5, 5, NULL);
    print_and_destroy(&g, &jdm);

    // Clean up
    igraph_destroy(&g);
    igraph_matrix_int_destroy(&jdm);
    igraph_vector_int_destroy(&weights);
    VERIFY_FINALLY_STACK();

    return 0;
}
