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

void call_and_print(igraph_t *graph, igraph_vector_bool_t *types) {
    igraph_matrix_t result;
    igraph_vector_t row_ids;
    igraph_vector_t col_ids;
    igraph_matrix_init(&result, 0, 0);
    igraph_vector_init(&row_ids, 0);
    igraph_vector_init(&col_ids, 0);
    IGRAPH_ASSERT(igraph_get_incidence(graph, types, &result, &row_ids, &col_ids) == IGRAPH_SUCCESS);
    printf("Incidence matrix:\n");
    print_matrix(&result);
    printf("Row ids:\n");
    print_vector(&row_ids);
    printf("Col ids:\n");
    print_vector(&col_ids);
    printf("\n");
    igraph_vector_destroy(&row_ids);
    igraph_vector_destroy(&col_ids);
    igraph_matrix_destroy(&result);
}


int main() {
    igraph_t g_0, g_1, g_mu, g_mun;
    igraph_vector_bool_t t_0, t_1, t_mu;
    igraph_matrix_t result;

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_small(&g_mu, 6, 0, 0,1, 0,2, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);
    igraph_small(&g_mun, 6, 0, 0,1, 0,2, 0,3, 1,3, 2,0, 2,0, 2,3, 3,4, 3,4, -1);

    igraph_vector_bool_init(&t_0, 0);
    igraph_vector_bool_init_int(&t_1, 1, 1);
    igraph_vector_bool_init_int(&t_mu, 6, 0, 1, 1, 0, 1, 0);

    igraph_matrix_init(&result, 0, 0);

    printf("No vertices:\n");
    call_and_print(&g_0, &t_0);

    printf("One vertex:\n");
    call_and_print(&g_1, &t_1);

    printf("Disconnected graph with multiple edges:\n");
    call_and_print(&g_mu, &t_mu);

    printf("Checking non-bipartite graph.\n");
    call_and_print(&g_mun, &t_mu);

    VERIFY_FINALLY_STACK();

    printf("Checking wrong type vector size error handling.\n");
    CHECK_ERROR(igraph_get_incidence(&g_mu, &t_0, &result, NULL, NULL), IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_mu);
    igraph_destroy(&g_mun);

    igraph_vector_bool_destroy(&t_0);
    igraph_vector_bool_destroy(&t_1);
    igraph_vector_bool_destroy(&t_mu);

    igraph_matrix_destroy(&result);

    VERIFY_FINALLY_STACK();
    return 0;
}
