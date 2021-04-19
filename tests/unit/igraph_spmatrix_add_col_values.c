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
    igraph_spmatrix_t spmat;
    int i;
    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("0x0 matrix, trying to add nonexistent column\n");
    igraph_spmatrix_init(&spmat, 0, 0);
    IGRAPH_ASSERT(igraph_spmatrix_add_col_values(&spmat, /*to*/ 0, /*from*/ 0) == IGRAPH_EINVAL);
    igraph_spmatrix_destroy(&spmat);

    printf("\n1x1 matrix, adding a column to itself\n");
    igraph_spmatrix_init(&spmat, 1, 1);
    igraph_spmatrix_set(&spmat, 0, 0, 5);
    IGRAPH_ASSERT(igraph_spmatrix_add_col_values(&spmat, /*to*/ 0, /*from*/ 0) == IGRAPH_SUCCESS);
    print_spmatrix(&spmat);
    igraph_spmatrix_destroy(&spmat);

    printf("\n5x6 matrix\n");
    igraph_spmatrix_init(&spmat, 5, 6);
    for (i = 0; i < 30; i++) {
        igraph_spmatrix_set(&spmat, i/6, i%6, i);
    }
    IGRAPH_ASSERT(igraph_spmatrix_add_col_values(&spmat, /*to*/ 1, /*from*/ 2) == IGRAPH_SUCCESS);
    print_spmatrix(&spmat);
    igraph_spmatrix_destroy(&spmat);

    VERIFY_FINALLY_STACK();
    return 0;
}
