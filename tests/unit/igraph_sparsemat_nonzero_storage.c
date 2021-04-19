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
    int i, j;
    int size = 3;

    printf("0x0 matrix\n");
    igraph_sparsemat_init(&spmat, 0, 0, 0);
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_nonzero_storage(&spmat) == 0);
    IGRAPH_ASSERT(igraph_sparsemat_nonzero_storage(&spmat_comp) == 0);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat_comp);

    printf("3x3 compressed matrix with duplicate values that add up to zero.\n");
    igraph_sparsemat_init(&spmat, size, size, 7);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            igraph_sparsemat_entry(&spmat, i, j, 5);
            igraph_sparsemat_entry(&spmat, i, j, -5);

            /* This checks if there's two entries for every loop. */
            IGRAPH_ASSERT(igraph_sparsemat_nonzero_storage(&spmat) == (i * size + j + 1) * 2);
            igraph_sparsemat_compress(&spmat, &spmat_comp);
            IGRAPH_ASSERT(igraph_sparsemat_nonzero_storage(&spmat_comp) == (i * size + j + 1) * 2);

            igraph_sparsemat_destroy(&spmat_comp);
        }
    }
    printf("Adding one entry to work around some broken error handling.\n");
    igraph_sparsemat_entry(&spmat, 0, 0, 5);
    igraph_sparsemat_compress(&spmat, &spmat_comp);

    printf("Removing duplicates should leave us with one entry in each position.\n");
    igraph_sparsemat_dupl(&spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_nonzero_storage(&spmat_comp) == (size * size));

    printf("Removing all zeros should leave us with only one entry.\n");
    igraph_sparsemat_dropzeros(&spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_nonzero_storage(&spmat_comp) == 1);

    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat_comp);

    VERIFY_FINALLY_STACK();
    return 0;
}
