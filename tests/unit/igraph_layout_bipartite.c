/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_matrix_t result;
    igraph_vector_bool_t types;

    printf("No vertices:\n");
    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init(&types, 0);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ 1.0, /*maxiter*/ 100) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("1 vertex:\n");
    igraph_small(&g, 1, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 1, 0);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ 1.0, /*maxiter*/ 100) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("4 vertices, disconnected, not actually bipartite, with loops and multiple edges:\n");
    igraph_small(&g, 4, 0, 0,0, 0,0, 0,0, 1,2, 1,2, 1,3, 1,3, 2,3, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 4, 0, 1, 0, 1);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ 1.0, /*maxiter*/ 100) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("10 vertices bipartite graph:\n");
    igraph_small(&g, 10, 0, 0,5, 0,7, 1,6, 1,7, 1,8, 2,5, 3,8, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 10, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ 1.0, /*maxiter*/100) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("10 vertices bipartite graph, no iterations:\n");
    igraph_small(&g, 10, 0, 0,5, 0,7, 1,6, 1,7, 1,8, 2,5, 3,8, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 10, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ 1.0, /*maxiter*/0) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("4 vertices with -10 true values for types:\n");
    igraph_small(&g, 4, 0, 0,1, 1,2, 2,3, 3,0, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 4, 0, -10, 0, -10);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ 1.0, /*maxiter*/ 100) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("4 vertices, negative vgaps:\n");
    igraph_small(&g, 4, 0, 0,1, 1,2, 2,3, 3,0, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 4, 0, 1, 0, 1);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ 1.0, /*vgap*/ -1.0, /*maxiter*/ 100) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("4 vertices, negative hgaps, emits error.\n");
    igraph_small(&g, 4, 0, 0,1, 1,2, 2,3, 3,0, -1);
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_bool_init_int(&types, 4, 0, 1, 0, 1);
    IGRAPH_ASSERT(igraph_layout_bipartite(&g, &types, &result, /*hgap*/ -1.0, /*vgap*/ 1.0, /*maxiter*/ 100) == IGRAPH_EINVAL);
    igraph_vector_bool_destroy(&types);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
