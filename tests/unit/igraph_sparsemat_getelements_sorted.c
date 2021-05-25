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

void print_res_i_j(igraph_vector_t *result, igraph_vector_int_t *i, igraph_vector_int_t *j) {
    print_vector(result);
    print_vector_int(i);
    print_vector_int(j);
}

void destroy_all(igraph_vector_t *result, igraph_vector_int_t *i, igraph_vector_int_t *j, igraph_sparsemat_t *spmat, igraph_sparsemat_t *spmat_comp) {
    igraph_vector_destroy(result);
    igraph_vector_int_destroy(i);
    igraph_vector_int_destroy(j);
    igraph_sparsemat_destroy(spmat);
    igraph_sparsemat_destroy(spmat_comp);
}

void init_all(igraph_vector_t *result, igraph_vector_int_t *i, igraph_vector_int_t *j, igraph_sparsemat_t *spmat) {
    igraph_vector_init(result, 0);
    igraph_vector_int_init(i, 0);
    igraph_vector_int_init(j, 0);
    igraph_sparsemat_init(spmat, 0, 0, 0);
}

int main() {
    igraph_sparsemat_t spmat;
    igraph_sparsemat_t spmat_comp;
    igraph_vector_t result;
    igraph_vector_int_t i, j;
    int k, l;
    int size = 3;

    printf("0x0 matrix\n");
    init_all(&result, &i, &j, &spmat);
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_getelements_sorted(&spmat, &i, &j, &result) == IGRAPH_SUCCESS);
    printf("triplet:\n");
    print_res_i_j(&result, &i, &j);
    IGRAPH_ASSERT(igraph_sparsemat_getelements_sorted(&spmat_comp, &i, &j, &result) == IGRAPH_SUCCESS);
    printf("compressed:\n");
    print_res_i_j(&result, &i, &j);
    destroy_all(&result, &i, &j, &spmat, &spmat_comp);

    printf("\n3x3 matrix\n");
    init_all(&result, &i, &j, &spmat);
    for (k = 0; k < size; k += 2) {
        for (l = 0; l < size; l ++) {
            igraph_sparsemat_entry(&spmat, k, l, 100);
            igraph_sparsemat_entry(&spmat, k, l, k * size + l);
        }
    }
    igraph_sparsemat_compress(&spmat, &spmat_comp);
    IGRAPH_ASSERT(igraph_sparsemat_getelements_sorted(&spmat, &i, &j, &result) == IGRAPH_SUCCESS);
    printf("triplet:\n");
    print_res_i_j(&result, &i, &j);
    IGRAPH_ASSERT(igraph_sparsemat_getelements_sorted(&spmat_comp, &i, &j, &result) == IGRAPH_SUCCESS);
    printf("compressed:\n");
    print_res_i_j(&result, &i, &j);

    destroy_all(&result, &i, &j, &spmat, &spmat_comp);

    VERIFY_FINALLY_STACK();
    return 0;
}
