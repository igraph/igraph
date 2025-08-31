/*
   igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_matrix_t mat;

    igraph_matrix_init(&mat, 3, 3);
    MATRIX(mat, 0, 0) =  1e-11;
    MATRIX(mat, 0, 1) = -1e-11;
    MATRIX(mat, 0, 2) =  1e-10;
    MATRIX(mat, 1, 0) = -1e-10;
    MATRIX(mat, 1, 1) = IGRAPH_INFINITY;
    MATRIX(mat, 1, 2) = -IGRAPH_INFINITY;
    MATRIX(mat, 2, 0) = IGRAPH_NAN;
    MATRIX(mat, 2, 1) = 0.0;
    MATRIX(mat, 2, 2) = 1.0;

    print_matrix(&mat);
    igraph_matrix_zapsmall(&mat, 0);
    print_matrix(&mat);

    CHECK_ERROR(igraph_matrix_zapsmall(&mat, -1), IGRAPH_EINVAL);

    igraph_matrix_destroy(&mat);

    /* TODO complex case */

    VERIFY_FINALLY_STACK();

    return 0;
}
