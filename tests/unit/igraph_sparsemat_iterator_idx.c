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
    igraph_sparsemat_iterator_t it;
    igraph_sparsemat_iterator_t it_comp;

    printf("0x0 matrix.\n");
    igraph_sparsemat_init(&spmat, 0, 0, /*nzmax*/0);
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    igraph_sparsemat_iterator_init(&it, &spmat);
    igraph_sparsemat_iterator_init(&it_comp, &spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_iterator_idx(&it) == 0);
    IGRAPH_ASSERT(igraph_sparsemat_iterator_idx(&it_comp) == 0);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat_comp);

    printf("3x3 matrix.\n");
    igraph_sparsemat_init(&spmat, 3, 3, /*nzmax*/0);
    igraph_sparsemat_entry(&spmat, 0, 0, 5);
    igraph_sparsemat_entry(&spmat, 3, 3, 6);
    igraph_sparsemat_entry(&spmat, 3, 3, 6);
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    igraph_sparsemat_iterator_init(&it, &spmat);
    igraph_sparsemat_iterator_init(&it_comp, &spmat_comp);
    igraph_sparsemat_iterator_next(&it);
    igraph_sparsemat_iterator_next(&it);
    igraph_sparsemat_iterator_next(&it_comp);
    IGRAPH_ASSERT(igraph_sparsemat_iterator_idx(&it) == 2);
    IGRAPH_ASSERT(igraph_sparsemat_iterator_idx(&it_comp) == 1);
    igraph_sparsemat_destroy(&spmat);
    igraph_sparsemat_destroy(&spmat_comp);

    VERIFY_FINALLY_STACK();
    return 0;
}
