/*
   IGraph library.
   Copyright (C) 2019-2022  The igraph development team <igraph@igraph.org>

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

#define SIMPLIFY_PRINT_DESTROY(name) \
    printf(name "\n"); \
    igraph_simplify_and_colorize(&graph, &res, &vcol, &ecol); \
    print_graph(&res); \
    print_vector_int(&vcol); \
    print_vector_int(&ecol); \
    printf("\n"); \
    igraph_destroy(&res); \
    igraph_destroy(&graph);

int main(void) {
    igraph_t graph, res;
    igraph_vector_int_t vcol, ecol;

    igraph_vector_int_init(&vcol, 0);
    igraph_vector_int_init(&ecol, 0);

    /* null graph */
    igraph_empty(&graph, 0, 0);
    SIMPLIFY_PRINT_DESTROY("K0");

    /* singleton graph */
    igraph_empty(&graph, 1, 0);
    SIMPLIFY_PRINT_DESTROY("K1");

    /* 4-cycle-graph */
    igraph_ring(&graph, 4, 0, 0, 1);
    SIMPLIFY_PRINT_DESTROY("C4");

    /* both multi-edges and self loops */
    igraph_small(&graph, 2, 0,
                 0, 1, 0, 1, 1, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Undirected graph 1");

    /* parallel edges specified with different vertex orderings */
    igraph_small(&graph, 3, 0,
                 0, 1, 1, 2, 2, 0, 2, 2, 2, 2, 2, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Undirected graph 2");

    /* directed version of the same as above */
    igraph_small(&graph, 3, 1,
                 0, 1, 1, 2, 2, 0, 2, 2, 2, 2, 2, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Directed graph 1");

    /* isolated vertices */
    igraph_small(&graph, 4, 1,
                 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, -1);
    SIMPLIFY_PRINT_DESTROY("Directed graph 2");

    igraph_vector_int_destroy(&vcol);
    igraph_vector_int_destroy(&ecol);

    VERIFY_FINALLY_STACK();

    return 0;
}
