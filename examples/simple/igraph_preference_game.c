/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <assert.h>

int main() {
    igraph_t g;
    igraph_vector_t tdist;
    igraph_matrix_t pmat;
    igraph_bool_t conn;
    igraph_vector_bool_t bs;
    int i;

    /* Symmetric preference game */
    igraph_vector_bool_init(&bs, 0);

    igraph_vector_init_real(&tdist, 3, 1.0, 1.0, 1.0);

    igraph_matrix_init(&pmat, 3, 3);
    for (i = 0; i < 3; i++) {
        MATRIX(pmat, i, i) = 0.2;
    }

    /* undirected, no loops */
    IGRAPH_CHECK(igraph_preference_game(&g, 1000, 3, &tdist, /*fixed_sizes=*/ 0,
                                        &pmat, 0, 0, 0));
    if (igraph_vcount(&g) != 1000) {
        return 18;
    }
    if (igraph_is_directed(&g)) {
        return 2;
    }
    igraph_is_connected(&g, &conn, IGRAPH_STRONG);
    if (conn) {
        return 3;
    }
    igraph_is_loop(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 4;
    }
    igraph_is_multiple(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 5;
    }
    igraph_destroy(&g);

    for (i = 0; i < 2; i++) {
        MATRIX(pmat, i, i + 1) = 0.1;
    }

    /* directed, no loops */
    IGRAPH_CHECK(igraph_preference_game(&g, 1000, 3, &tdist, /*fixed_sizes=*/0,
                                        &pmat, 0, 1, 0));
    if (igraph_vcount(&g) != 1000) {
        return 17;
    }
    if (!igraph_is_directed(&g)) {
        return 6;
    }
    igraph_is_loop(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 7;
    }
    igraph_is_multiple(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 8;
    }
    igraph_destroy(&g);

    /* undirected, loops */
    for (i = 0; i < 3; i++) {
        MATRIX(pmat, i, i) = 1.0;
    }
    IGRAPH_CHECK(igraph_preference_game(&g, 100, 3, &tdist, /*fixed_sizes=*/ 0,
                                        &pmat, 0, 0, 1));
    if (igraph_vcount(&g) != 100) {
        return 16;
    }
    if (igraph_ecount(&g) < 1395) {
        return 20;
    }
    if (igraph_is_directed(&g)) {
        return 9;
    }
    igraph_is_loop(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs) == 0) {
        return 10;
    }
    igraph_is_multiple(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 11;
    }
    igraph_destroy(&g);

    /* directed, loops */
    IGRAPH_CHECK(igraph_preference_game(&g, 100, 3, &tdist, /*fixed_sizes=*/ 0,
                                        &pmat, 0, 1, 1));
    if (igraph_vcount(&g) != 100) {
        return 15;
    }
    if (igraph_ecount(&g) < 2700) {
        return 19;
    }
    if (!igraph_is_directed(&g)) {
        return 12;
    }
    igraph_is_loop(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs) == 0) {
        return 13;
    }
    igraph_is_multiple(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 14;
    }
    igraph_destroy(&g);

    /* Asymmetric preference game */

    /* directed, no loops */
    igraph_matrix_resize(&pmat, 2, 2);
    MATRIX(pmat, 0, 0) = 1;
    MATRIX(pmat, 0, 1) = 1;
    MATRIX(pmat, 1, 0) = 1;
    MATRIX(pmat, 1, 1) = 1;
    IGRAPH_CHECK(igraph_asymmetric_preference_game(&g, 100, 2, 0, &pmat, 0, 0, 0));
    if (igraph_vcount(&g) != 100) {
        return 21;
    }
    if (igraph_ecount(&g) != 9900) {
        return 22;
    }
    if (!igraph_is_directed(&g)) {
        return 23;
    }
    igraph_is_loop(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 24;
    }
    igraph_is_multiple(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 25;
    }
    igraph_destroy(&g);

    /* directed, loops */
    igraph_matrix_resize(&pmat, 2, 2);
    MATRIX(pmat, 0, 0) = 1;
    MATRIX(pmat, 0, 1) = 1;
    MATRIX(pmat, 1, 0) = 1;
    MATRIX(pmat, 1, 1) = 1;
    IGRAPH_CHECK(igraph_asymmetric_preference_game(&g, 100, 2, 0, &pmat, 0, 0, 1));
    if (igraph_vcount(&g) != 100) {
        return 26;
    }
    if (igraph_ecount(&g) != 10000) {
        return 27;
    }
    if (!igraph_is_directed(&g)) {
        return 28;
    }
    igraph_is_loop(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs) != 100) {
        return 29;
    }
    igraph_is_multiple(&g, &bs, igraph_ess_all(IGRAPH_EDGEORDER_ID));
    if (igraph_vector_bool_sum(&bs)) {
        return 30;
    }
    igraph_destroy(&g);

    igraph_vector_destroy(&tdist);
    igraph_matrix_destroy(&pmat);
    igraph_vector_bool_destroy(&bs);

    assert(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;
}
