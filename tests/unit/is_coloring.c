/*
   igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

void test_vertex_coloring(void) {
    igraph_t graph;
    igraph_vector_int_t types;
    igraph_bool_t res;

    printf("Testing igraph_is_vertex_coloring()\n");

    /* Empty graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&types, 0);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Single vertex */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_init_int(&types, 1, 0);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Two connected vertices with different colors */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_vector_int_init_int(&types, 2, 0, 1);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);

    /* Two connected vertices with same colors (invalid) */
    igraph_vector_int_init_int(&types, 2, 0, 0);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(!res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Triangle with valid 3-coloring */
    igraph_small(&graph, 3, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_vector_int_init_int(&types, 3, 0, 1, 2);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);

    /* Triangle with invalid coloring (two adjacent vertices same color) */
    igraph_vector_int_init_int(&types, 3, 0, 0, 1);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(!res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Self-loop should not affect coloring validity */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 0, 0, 0, 1, -1);
    igraph_vector_int_init_int(&types, 2, 0, 1);
    igraph_is_vertex_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    printf("igraph_is_vertex_coloring() tests passed\n");
}

void test_bipartite_coloring(void) {
    igraph_t graph;
    igraph_vector_bool_t types;
    igraph_bool_t res;
    igraph_neimode_t mode;

    printf("Testing igraph_is_bipartite_coloring()\n");

    /* Empty graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_vector_bool_init(&types, 0);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(res);
    IGRAPH_ASSERT(mode == IGRAPH_ALL);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Single vertex */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_vector_bool_init_int(&types, 1, 0);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(res);
    IGRAPH_ASSERT(mode == IGRAPH_ALL);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Valid bipartite undirected graph */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 2, 0, 3, 1, 2, 1, 3, -1);
    igraph_vector_bool_init_int(&types, 4, 0, 0, 1, 1);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(res);
    IGRAPH_ASSERT(mode == IGRAPH_ALL);  /* Undirected graph */
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Invalid bipartite coloring (two adjacent vertices same type) */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 1, 3, 2, 3, -1);
    igraph_vector_bool_init_int(&types, 4, 0, 0, 1, 1);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(!res);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Directed bipartite graph - all edges from false to true */
    igraph_small(&graph, 4, IGRAPH_DIRECTED, 0, 2, 0, 3, 1, 2, 1, 3, -1);
    igraph_vector_bool_init_int(&types, 4, 0, 0, 1, 1);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(res);
    IGRAPH_ASSERT(mode == IGRAPH_OUT);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Directed bipartite graph - all edges from true to false */
    igraph_small(&graph, 4, IGRAPH_DIRECTED, 2, 0, 3, 0, 2, 1, 3, 1, -1);
    igraph_vector_bool_init_int(&types, 4, 0, 0, 1, 1);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(res);
    IGRAPH_ASSERT(mode == IGRAPH_IN);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Directed bipartite graph - edges in both directions */
    igraph_small(&graph, 4, IGRAPH_DIRECTED, 0, 2, 2, 1, 1, 3, 3, 0, -1);
    igraph_vector_bool_init_int(&types, 4, 0, 0, 1, 1);
    igraph_is_bipartite_coloring(&graph, &types, &res, &mode);
    IGRAPH_ASSERT(res);
    IGRAPH_ASSERT(mode == IGRAPH_ALL);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    /* Self-loop should not affect bipartite coloring validity */
    igraph_small(&graph, 3, IGRAPH_UNDIRECTED, 0, 0, 0, 2, 1, 2, -1);
    igraph_vector_bool_init_int(&types, 3, 0, 0, 1);
    igraph_is_bipartite_coloring(&graph, &types, &res, NULL);
    IGRAPH_ASSERT(res);
    igraph_vector_bool_destroy(&types);
    igraph_destroy(&graph);

    printf("igraph_is_bipartite_coloring() tests passed\n");
}

void test_edge_coloring(void) {
    igraph_t graph;
    igraph_vector_int_t types;
    igraph_bool_t res;

    printf("Testing igraph_is_edge_coloring()\n");

    /* Empty graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&types, 0);
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Single vertex, no edges */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_vector_int_init(&types, 0);
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Two vertices, one edge */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 0, 1, -1);
    igraph_vector_int_init_int(&types, 1, 0);
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Triangle with valid edge coloring */
    igraph_small(&graph, 3, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_vector_int_init_int(&types, 3, 0, 1, 2);  /* Three different colors */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);

    /* Triangle with invalid edge coloring (two incident edges same color) */
    igraph_vector_int_init_int(&types, 3, 0, 1, 0);  /* Edges 0 and 2 both color 0, both incident to vertex 0 */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(!res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Star graph with valid edge coloring */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 0, 3, -1);
    igraph_vector_int_init_int(&types, 3, 0, 1, 2);  /* All edges incident to vertex 0 have different colors */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);

    /* Star graph with invalid edge coloring */
    igraph_vector_int_init_int(&types, 3, 0, 0, 1);  /* Two edges incident to vertex 0 have same color */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(!res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Path graph with valid edge coloring */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, -1);
    igraph_vector_int_init_int(&types, 3, 0, 1, 0);  /* Alternating colors */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    /* Graph with self-loop: 1 - 1 - 2 (self-loop and regular edge) */
    /* Self-loops are not considered adjacent to themselves */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 1, 1, 1, 0, -1);  /* vertex 1 self-loop, edge 1-0 */
    igraph_vector_int_init_int(&types, 2, 0, 1);  /* Different colors - should be valid */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(res);
    igraph_vector_int_destroy(&types);

    /* Same colors for self-loop and adjacent edge - should be invalid */
    igraph_vector_int_init_int(&types, 2, 0, 0);  /* Same colors - should be invalid */
    igraph_is_edge_coloring(&graph, &types, &res);
    IGRAPH_ASSERT(!res);
    igraph_vector_int_destroy(&types);
    igraph_destroy(&graph);

    printf("igraph_is_edge_coloring() tests passed\n");
}

void test_error_conditions(void) {
    igraph_t graph;
    igraph_vector_int_t int_types;
    igraph_vector_bool_t bool_types;
    igraph_bool_t res;
    igraph_neimode_t mode;

    printf("Testing error conditions\n");

    /* Create a simple graph */
    igraph_small(&graph, 3, IGRAPH_UNDIRECTED, 0, 1, 1, 2, -1);

    /* Wrong size vector for vertex coloring */
    igraph_vector_int_init_int(&int_types, 2, 0, 1);  /* Size 2, but graph has 3 vertices */
    CHECK_ERROR(igraph_is_vertex_coloring(&graph, &int_types, &res), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&int_types);

    /* Wrong size vector for bipartite coloring */
    igraph_vector_bool_init_int(&bool_types, 2, 0, 1);  /* Size 2, but graph has 3 vertices */
    CHECK_ERROR(igraph_is_bipartite_coloring(&graph, &bool_types, &res, &mode), IGRAPH_EINVAL);
    igraph_vector_bool_destroy(&bool_types);

    igraph_destroy(&graph);

    /* Wrong size vector for edge coloring */
    igraph_small(&graph, 2, IGRAPH_UNDIRECTED, 0, 1, -1);  /* Graph with 1 edge */
    igraph_vector_int_init_int(&int_types, 2, 0, 1);  /* Size 2, but graph has 1 edge */
    CHECK_ERROR(igraph_is_edge_coloring(&graph, &int_types, &res), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&int_types);
    igraph_destroy(&graph);

    printf("Error condition tests passed\n");
}

int main(void) {
    test_vertex_coloring();
    test_bipartite_coloring();
    test_edge_coloring();
    test_error_conditions();

    VERIFY_FINALLY_STACK();

    return 0;
}
