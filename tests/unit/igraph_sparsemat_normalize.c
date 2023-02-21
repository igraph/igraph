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
#include "test_utilities.h"

void create_test_matrix(
    igraph_sparsemat_t* spmat,
    igraph_sparsemat_t* spmat_comp,
    igraph_bool_t use_zero_col_row
) {
    igraph_sparsemat_init(spmat, 3, 3, /*nzmax*/7);
    igraph_sparsemat_entry(spmat, 0, 0, 6);
    igraph_sparsemat_entry(spmat, 2, 2, 7);

    if (!use_zero_col_row) {
        igraph_sparsemat_entry(spmat, 0, 1, 2);
        igraph_sparsemat_entry(spmat, 1, 0, 4);
        igraph_sparsemat_entry(spmat, 1, 1, 6);
        igraph_sparsemat_entry(spmat, 2, 1, 2);
        igraph_sparsemat_entry(spmat, 1, 2, -14);
    }

    igraph_sparsemat_compress(spmat, spmat_comp);
}

int main(void) {
    igraph_sparsemat_t spmat;
    igraph_sparsemat_t spmat_comp;

    printf("0x0 matrix.\n");
    igraph_sparsemat_init(&spmat, 0, 0, /*nzmax*/0);
    igraph_sparsemat_compress(&spmat, &spmat_comp);

    IGRAPH_ASSERT(igraph_sparsemat_normalize_rows(&spmat_comp, /* allow_zero = */ 1) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_cols(&spmat_comp, /* allow_zero = */ 1) == IGRAPH_SUCCESS);

    /* Normalization should succeed for an empty matrix even if allow_zero = 0
     * because technically there is no 0/0 division anywhere */
    IGRAPH_ASSERT(igraph_sparsemat_normalize_rows(&spmat_comp, /* allow_zero = */ 0) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_cols(&spmat_comp, /* allow_zero = */ 0) == IGRAPH_SUCCESS);

    igraph_sparsemat_destroy(&spmat_comp);
    igraph_sparsemat_destroy(&spmat);

    printf("\n");
    printf("3x3 matrix without zero rows and columns, row-wise normalization.\n");
    create_test_matrix(&spmat, &spmat_comp, /* use_zero_col_row = */ 0);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_rows(&spmat_comp, /* allow_zero = */ 0) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    igraph_sparsemat_destroy(&spmat_comp);
    igraph_sparsemat_destroy(&spmat);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("\n");
    printf("3x3 matrix without zero rows and columns, column-wise normalization.\n");
    create_test_matrix(&spmat, &spmat_comp, /* use_zero_col_row = */ 0);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_cols(&spmat_comp, /* allow_zero = */ 0) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    igraph_sparsemat_destroy(&spmat_comp);
    igraph_sparsemat_destroy(&spmat);

    printf("\n");
    printf("3x3 matrix with zero rows and columns, row-wise normalization.\n");
    create_test_matrix(&spmat, &spmat_comp, /* use_zero_col_row = */ 1);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_rows(&spmat_comp, /* allow_zero = */ 0) == IGRAPH_EINVAL);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_rows(&spmat_comp, /* allow_zero = */ 1) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    igraph_sparsemat_destroy(&spmat_comp);
    igraph_sparsemat_destroy(&spmat);

    printf("\n");
    printf("3x3 matrix with zero rows and columns, column-wise normalization.\n");
    create_test_matrix(&spmat, &spmat_comp, /* use_zero_col_row = */ 1);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_cols(&spmat_comp, /* allow_zero = */ 0) == IGRAPH_EINVAL);
    IGRAPH_ASSERT(igraph_sparsemat_normalize_cols(&spmat_comp, /* allow_zero = */ 1) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    igraph_sparsemat_destroy(&spmat_comp);
    igraph_sparsemat_destroy(&spmat);

    printf("\n");
    printf("uncompressed matrix.\n");
    igraph_sparsemat_init(&spmat, 0, 0, /*nzmax*/0);
    IGRAPH_ASSERT(igraph_sparsemat_droptol(&spmat, 10) == IGRAPH_EINVAL);
    igraph_sparsemat_destroy(&spmat);

    VERIFY_FINALLY_STACK();

    return 0;
}
