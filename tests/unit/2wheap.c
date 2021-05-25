/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

#include "core/indheap.h"

#include "test_utilities.inc"

int main() {

    igraph_vector_t elems;
    igraph_2wheap_t Q;
    long int i;
    igraph_real_t prev = IGRAPH_INFINITY;

    srand(42); /* make tests deterministic */

    igraph_vector_init(&elems, 100);
    for (i = 0; i < igraph_vector_size(&elems); i++) {
        VECTOR(elems)[i] = rand() / (double)RAND_MAX;
    }

    igraph_2wheap_init(&Q, igraph_vector_size(&elems));
    for (i = 0; i < igraph_vector_size(&elems); i++) {
        igraph_2wheap_push_with_index(&Q, i, VECTOR(elems)[i]);
    }

    /*****/

    for (i = 0; i < igraph_vector_size(&elems); i++) {
        if (VECTOR(elems)[i] != igraph_2wheap_get(&Q, i)) {
            return 1;
        }
    }

    /*****/

    for (i = 0; i < igraph_vector_size(&elems); i++) {
        long int j;
        igraph_real_t tmp = igraph_2wheap_max(&Q);
        if (tmp > prev) {
            return 2;
        }
        if (tmp != igraph_2wheap_delete_max_index(&Q, &j)) {
            return 3;
        }
        if (VECTOR(elems)[j] != tmp) {
            return 4;
        }
        prev = tmp;
    }

    /*****/

    for (i = 0; i < igraph_vector_size(&elems); i++) {
        igraph_2wheap_push_with_index(&Q, i, VECTOR(elems)[i]);
    }
    if (igraph_2wheap_size(&Q) != igraph_vector_size(&elems)) {
        return 5;
    }
    for (i = 0; i < igraph_vector_size(&elems); i++) {
        VECTOR(elems)[i] = rand() / (double)RAND_MAX;
        igraph_2wheap_modify(&Q, i, VECTOR(elems)[i]);
    }
    for (i = 0; i < igraph_vector_size(&elems); i++) {
        if (VECTOR(elems)[i] != igraph_2wheap_get(&Q, i)) {
            return 6;
        }
    }
    prev = IGRAPH_INFINITY;
    for (i = 0; i < igraph_vector_size(&elems); i++) {
        long int j;
        igraph_real_t tmp = igraph_2wheap_max(&Q);
        if (tmp > prev) {
            return 7;
        }
        if (tmp != igraph_2wheap_delete_max_index(&Q, &j)) {
            return 8;
        }
        if (VECTOR(elems)[j] != tmp) {
            return 9;
        }
        prev = tmp;
    }
    if (!igraph_2wheap_empty(&Q)) {
        return 10;
    }
    if (igraph_2wheap_size(&Q) != 0) {
        return 11;
    }

    igraph_2wheap_destroy(&Q);
    igraph_vector_destroy(&elems);

    /* Hand-made example */

#define MAX       do { igraph_2wheap_delete_max(&Q); igraph_2wheap_check(&Q); } while (0)
#define PUSH(i,e) do { igraph_2wheap_push_with_index(&Q, (i), -(e)); igraph_2wheap_check(&Q); } while (0);
#define MOD(i, e) do { igraph_2wheap_modify(&Q, (i), -(e)); igraph_2wheap_check(&Q); } while (0)

    igraph_2wheap_init(&Q, 21);
    /* 0.00 [ 4] */ PUSH(4, 0);
    /* MAX       */ MAX;
    /* 0.63 [11] */ PUSH(11, 0.63);
    /* 0.05 [15] */ PUSH(15, 0.05);
    /* MAX       */ MAX;
    /* 0.4  [12] */ PUSH(12, 0.4);
    /* 0.4  [13] */ PUSH(13, 0.4);
    /* 0.12 [16] */ PUSH(16, 0.12);
    /* MAX       */ MAX;
    /* 1.1  [ 0] */ PUSH(0, 1.1);
    /* 1.1  [14] */ PUSH(14, 1.1);
    /* MAX       */ MAX;
    /* [11]/0.44 */ MOD(11, 0.44);
    /* MAX       */ MAX;
    /* MAX       */ MAX;
    /* 1.1  [20] */ PUSH(20, 1.1);
    /* MAX       */ MAX;
    /* 1.3  [ 7] */ PUSH(7, 1.3);
    /* 1.7  [ 9] */ PUSH(9, 1.7);
    /* MAX       */ MAX;
    /* 1.6  [19] */ PUSH(19, 1.6);
    /* MAX       */ MAX;
    /* 2.1  [17] */ PUSH(17, 2.1);
    /* 1.3  [18] */ PUSH(18, 1.3);
    /* MAX       */ MAX;
    /* 2.3  [ 1] */ PUSH(1, 2.3);
    /* 2.2  [ 5] */ PUSH(5, 2.2);
    /* 2.3  [10] */ PUSH(10, 2.3);
    /* MAX       */ MAX;
    /* [17]/1.5  */ MOD(17, 1.5);
    /* MAX       */ MAX;
    /* 1.8  [ 6] */ PUSH(6, 1.8);
    /* MAX       */ MAX;
    /* 1.3  [ 3] */ PUSH(3, 1.3);
    /* [ 6]/1.3  */ MOD(6, 1.3);
    /* MAX       */ MAX;
    /* 1.6  [ 8] */ PUSH(8, 1.6);
    /* MAX       */ MAX;

    igraph_2wheap_destroy(&Q);

    VERIFY_FINALLY_STACK();

    return 0;
}
