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
    igraph_sparsemat_t spmat_comp;

    printf("0x0 matrix.\n");
    igraph_sparsemat_init(&spmat, 0, 0, /*nzmax*/0);
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_droptol(&spmat_comp, 5) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat_comp);

    printf("3x3 matrix.\n");
    igraph_sparsemat_init(&spmat, 3, 3, /*nzmax*/7);
    igraph_sparsemat_entry(&spmat, 0, 0, 5);
    igraph_sparsemat_entry(&spmat, 1, 1, 6);
    igraph_sparsemat_entry(&spmat, 2, 2, 7);
    igraph_sparsemat_entry(&spmat, 3, 0, 1);
    igraph_sparsemat_entry(&spmat, 0, 3, 2);
    igraph_sparsemat_entry(&spmat, 2, 1, 3);
    igraph_sparsemat_entry(&spmat, 1, 2, -14);
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    printf("Remove values within distance 5 from zero:\n");
    IGRAPH_ASSERT(igraph_sparsemat_droptol(&spmat_comp, 5) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    printf("Remove values within distance 20 from zero:\n");
    IGRAPH_ASSERT(igraph_sparsemat_droptol(&spmat_comp, 20) == IGRAPH_SUCCESS);
    igraph_sparsemat_print(&spmat_comp, stdout);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat_comp);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("uncompressed matrix.\n");
    igraph_sparsemat_init(&spmat, 0, 0, /*nzmax*/0);
    IGRAPH_ASSERT(igraph_sparsemat_droptol(&spmat, 10) == IGRAPH_EINVAL);
    igraph_sparsemat_destroy(&spmat);

    VERIFY_FINALLY_STACK();
    return 0;
}
