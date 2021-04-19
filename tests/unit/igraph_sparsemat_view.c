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
    igraph_sparsemat_t spmat;
    igraph_matrix_t mat;
    int p[] = {0, 1, 3};
    int i[]  = {1, 0, 2};
    double x[] = {1, 5, 2};

    printf("Empty sparsemat.\n");
    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_view(&spmat, /*nzmax*/ 0, /*m*/ 0, /*n*/ 0, /*p*/ NULL, /*i*/ NULL, /*x*/ NULL, /*nz*/ 0);
    igraph_sparsemat_as_matrix(&mat, &spmat);
    print_matrix(&mat);
    igraph_free(spmat.cs);
    igraph_matrix_destroy(&mat);

    printf("3x2 sparsemat:\n");
    igraph_matrix_init(&mat, 0, 0);
    igraph_sparsemat_view(&spmat, /*nzmax*/ 3, /*m*/ 3, /*n*/ 2, p, i, x, /*nz*/ -1);
    igraph_sparsemat_as_matrix(&mat, &spmat);
    print_matrix(&mat);
    igraph_free(spmat.cs);
    igraph_matrix_destroy(&mat);

    VERIFY_FINALLY_STACK();
    return 0;
}
