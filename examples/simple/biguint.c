/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include "igraph_types_internal.h"
#include "bigint.h"

#include <limits.h>

int main() {

    igraph_biguint_t A, B, C, D, E, zero, one;

    igraph_biguint_init(&A);
    igraph_biguint_init(&B);
    igraph_biguint_init(&C);
    igraph_biguint_init(&D);
    igraph_biguint_init(&E);
    igraph_biguint_init(&zero);
    igraph_biguint_init(&one);

    /* set & add & sub */
    igraph_biguint_set_limb(&one, 1);
    igraph_biguint_set_limb(&A, UINT_MAX);
    igraph_biguint_set_limb(&B, UINT_MAX);
    igraph_biguint_add(&A, &A, &B);                      /* A <- A + B */

    igraph_biguint_print(&B);
    putchar('\n');
    igraph_biguint_print(&A);
    putchar('\n');

    igraph_biguint_sub(&A, &A, &B);                      /* A <- A - B */
    if (!igraph_biguint_equal(&A, &B)) {
        return 1;
    }

    /* inc & dec */
    igraph_biguint_inc(&A, &A);                          /* A <- A + 1 */
    igraph_biguint_dec(&A, &A);                          /* A <- A - 1 */
    if (!igraph_biguint_equal(&A, &B)) {
        return 2;
    }

    /* mul & div */
    igraph_biguint_mul(&C, &A, &B);                      /* C <- A * B */
    igraph_biguint_div(&E, &D, &C, &B);                  /* E <- C / B */
    /* D <- C % B */
    if (!igraph_biguint_equal(&E, &A)) {
        return 3;
    }
    if (!igraph_biguint_equal(&D, &zero)) {
        return 4;
    }

    igraph_biguint_mul(&C, &A, &A);                      /* C <- A * A */
    igraph_biguint_mul(&D, &C, &A);              /* C <- C * A */
    igraph_biguint_mul(&C, &D, &A);              /* C <- C * A */

    igraph_biguint_div(&C, &D, &C, &A);                  /* C <- C / A */
    igraph_biguint_div(&C, &D, &C, &A);              /* C <- C / A */
    igraph_biguint_div(&C, &D, &C, &A);              /* C <- C / A */
    igraph_biguint_div(&C, &D, &C, &A);              /* C <- C / A */
    if (!igraph_biguint_equal(&C, &one)) {
        return 5;
    }

    igraph_biguint_destroy(&A);
    igraph_biguint_destroy(&B);
    igraph_biguint_destroy(&C);
    igraph_biguint_destroy(&D);
    igraph_biguint_destroy(&E);
    igraph_biguint_destroy(&zero);
    igraph_biguint_destroy(&one);

    return 0;
}
