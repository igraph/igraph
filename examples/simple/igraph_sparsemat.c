/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int main() {

    igraph_sparsemat_t A, B, C, D;
    igraph_t G, H;
    igraph_vector_t vect;
    long int i;

    /* Create, compress, destroy */
    igraph_sparsemat_init(&A, 100, 20, 50);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_destroy(&B);
    igraph_sparsemat_destroy(&A);

    /* Convert a ring graph to a matrix, print it, compress, print again */
#define VC 10
    igraph_ring(&G, VC, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_get_sparsemat(&G, &A);
    igraph_destroy(&G);

    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_print(&A, stdout);
    igraph_sparsemat_print(&B, stdout);

    /* Basic query, nrow, ncol, type, is_triplet, is_cc */
    if (igraph_sparsemat_nrow(&A) != VC ||
        igraph_sparsemat_ncol(&A) != VC ||
        igraph_sparsemat_nrow(&B) != VC ||
        igraph_sparsemat_ncol(&B) != VC) {
        return 1;
    }
    if (!igraph_sparsemat_is_triplet(&A)) {
        return 2;
    }
    if (!igraph_sparsemat_is_cc(&B))      {
        return 3;
    }
    if (igraph_sparsemat_type(&A) != IGRAPH_SPARSEMAT_TRIPLET) {
        return 4;
    }
    if (igraph_sparsemat_type(&B) != IGRAPH_SPARSEMAT_CC)      {
        return 5;
    }

    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);
#undef VC

    printf("------------------------\n");

    /* Create unit matrices */
    igraph_sparsemat_eye(&A, /*n=*/ 5, /*nzmax=*/ 5, /*value=*/ 1.0,
                         /*compress=*/ 0);
    igraph_sparsemat_eye(&B, /*n=*/ 5, /*nzmax=*/ 5, /*value=*/ 1.0,
                         /*compress=*/ 1);
    igraph_sparsemat_print(&A, stdout);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);

    printf("------------------------\n");

    /* Create diagonal matrices */
    igraph_vector_init(&vect, 5);
    for (i = 0; i < 5; i++) {
        VECTOR(vect)[i] = i;
    }
    igraph_sparsemat_diag(&A, /*nzmax=*/ 5, /*values=*/ &vect, /*compress=*/ 0);
    igraph_sparsemat_diag(&B, /*nzmax=*/ 5, /*values=*/ &vect, /*compress=*/ 1);
    igraph_vector_destroy(&vect);
    igraph_sparsemat_print(&A, stdout);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);

    printf("------------------------\n");

    /* Transpose matrices */
    igraph_tree(&G, 10, /*children=*/ 2, IGRAPH_TREE_OUT);
    igraph_get_sparsemat(&G, &A);
    igraph_destroy(&G);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_transpose(&B, &C, /*values=*/ 1);
    igraph_sparsemat_print(&C, stdout);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);
    igraph_sparsemat_destroy(&C);

    printf("------------------------\n");

    /* Add duplicate elements */
    igraph_sparsemat_init(&A, 10, 10, /*nzmax=*/ 20);
    for (i = 1; i < 10; i++) {
        igraph_sparsemat_entry(&A, 0, i, 1.0);
    }
    for (i = 1; i < 10; i++) {
        igraph_sparsemat_entry(&A, 0, i, 1.0);
    }
    igraph_sparsemat_print(&A, stdout);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_dupl(&B);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);

    printf("------------------------\n");

    /* Drop zero elements */
    igraph_sparsemat_init(&A, 10, 10, /*nzmax=*/ 20);
    igraph_sparsemat_entry(&A, 7, 3, 0.0);
    for (i = 1; i < 10; i++) {
        igraph_sparsemat_entry(&A, 0, i, 1.0);
        igraph_sparsemat_entry(&A, 0, i, 0.0);
    }
    igraph_sparsemat_entry(&A, 0, 0, 0.0);
    igraph_sparsemat_print(&A, stdout);
    igraph_sparsemat_compress(&A, &B);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_dropzeros(&B);
    igraph_sparsemat_print(&B, stdout);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);

    printf("------------------------\n");

    /* Add two matrices */

    igraph_star(&G, 10, IGRAPH_STAR_OUT, /*center=*/ 0);
    igraph_ring(&H, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
    igraph_get_sparsemat(&G, &A);
    igraph_get_sparsemat(&H, &B);
    igraph_destroy(&G);
    igraph_destroy(&H);
    igraph_sparsemat_compress(&A, &C);
    igraph_sparsemat_compress(&B, &D);
    igraph_sparsemat_destroy(&A);
    igraph_sparsemat_destroy(&B);
    igraph_sparsemat_add(&C, &D, /*alpha=*/ 1.0, /*beta=*/ 2.0, &A);
    igraph_sparsemat_destroy(&C);
    igraph_sparsemat_destroy(&D);
    igraph_sparsemat_print(&A, stdout);
    igraph_sparsemat_destroy(&A);

    return 0;
}
