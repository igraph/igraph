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
    igraph_vector_int_t result;
    igraph_matrix_t m_pdc, m_0, m_m34, m_m43;
    int pdc[] = {9, 2, 7, 8,
                 6, 4, 3, 7,
                 5, 8, 1, 8,
                 7, 6, 9, 4};
    int m34[] = {3, 3, 2, 3,
                 2, 3, 3, 3,
                 3, 2, 3, 3};
    int m43[] = {3, 3, 2,
                 2, 3, 3,
                 3, 2, 3,
                 2, 3, 3};
    igraph_vector_int_init(&result, 0);
    matrix_init_int_row_major(&m_pdc, 4, 4, pdc);
    matrix_init_int_row_major(&m_m34, 3, 4, m34);
    matrix_init_int_row_major(&m_m43, 4, 3, m43);
    igraph_matrix_init(&m_0, 0, 0);

    printf("4 tasks, 4 agents:\n");
    igraph_solve_lsap(&m_pdc, 4, &result);
    print_vector_int(&result);

    printf("\n0 tasks, 0 agents:\n");
    igraph_solve_lsap(&m_0, 0, &result);
    print_vector_int(&result);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("\n4 tasks, 3 agents, n = 4.\n");
    IGRAPH_ASSERT(igraph_solve_lsap(&m_m34, 4, &result) == IGRAPH_EINVAL);

    printf("\n3 tasks, 4 agents, n = 4.\n");
    IGRAPH_ASSERT(igraph_solve_lsap(&m_m43, 4, &result) == IGRAPH_EINVAL);

    igraph_vector_int_destroy(&result);
    igraph_matrix_destroy(&m_pdc);
    igraph_matrix_destroy(&m_0);
    igraph_matrix_destroy(&m_m34);
    igraph_matrix_destroy(&m_m43);

    VERIFY_FINALLY_STACK();
    return 0;
}
