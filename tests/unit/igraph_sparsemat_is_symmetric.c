/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "test_utilities.inc"

#define DIM 10

#define INT(a) (igraph_rng_get_integer(igraph_rng_default(), 0, (a)))

int main() {
    int runs = 100;
    const int noelements = 20;
    igraph_sparsemat_t A;
    int i;

    igraph_rng_seed(igraph_rng_default(), 42);

    for (; runs > 0; runs--) {

        igraph_sparsemat_init(&A, DIM, DIM, noelements * 2);
        for (i = 0; i < noelements; i++) {
            int row = INT(DIM - 1);
            int col = INT(DIM - 1);
            int val = INT(100);
            igraph_sparsemat_entry(&A, row, col, val);
            igraph_sparsemat_entry(&A, col, row, val);
        }
        if (!igraph_sparsemat_is_symmetric(&A)) {
            return 1;
        }
        igraph_sparsemat_destroy(&A);

        igraph_sparsemat_init(&A, DIM, DIM, noelements);
        for (i = 0; i < noelements; i++) {
            igraph_sparsemat_entry(&A, INT(DIM - 1), INT(DIM - 1), INT(100));
        }
        if (igraph_sparsemat_is_symmetric(&A)) {
            return 2;
        }
        igraph_sparsemat_destroy(&A);

    }

    VERIFY_FINALLY_STACK();

    return 0;
}
