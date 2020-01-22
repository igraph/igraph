/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA 02139, USA

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
    igraph_vector_t outseq, inseq;
    igraph_bool_t result;

    /***** Testing igraph_is_degree_sequence *****/

    /* Valid undirected degree sequence */
    igraph_vector_init_int_end(&outseq, -1, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    igraph_is_degree_sequence(&outseq, 0, &result);
    if (!result) {
        return 1;
    }
    igraph_vector_destroy(&outseq);

    /* Undirected degree sequence with negative degree */
    igraph_vector_init_int_end(&outseq, -1, 3, -2, 3, 3, 3, 3, 3, 3, -1);
    igraph_is_degree_sequence(&outseq, 0, &result);
    if (result) {
        return 2;
    }
    igraph_vector_destroy(&outseq);

    /* Undirected degree sequence with uneven sum */
    igraph_vector_init_int_end(&outseq, -1, 3, 3, 3, 3, 3, 3, 3, -1);
    igraph_is_degree_sequence(&outseq, 0, &result);
    if (result) {
        return 3;
    }
    igraph_vector_destroy(&outseq);

    /* Valid directed degree sequences */
    igraph_vector_init_int_end(&outseq, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_init_int_end(&inseq, -1, 0, 3, 1, 3, 2, 4, 4, 1, 3, 1, -1);
    igraph_is_degree_sequence(&outseq, &inseq, &result);
    if (!result) {
        return 4;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Directed degree sequence with negative degree */
    igraph_vector_init_int_end(&outseq, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_init_int_end(&inseq, -1, 0, 3, 1, -7, 2, 4, 4, 1, 3, 1, -1);
    igraph_is_degree_sequence(&outseq, &inseq, &result);
    if (result) {
        return 5;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Directed degree sequence with different lengths */
    igraph_vector_init_int_end(&outseq, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_init_int_end(&inseq, -1, 0, 3, 1, 2, 4, 4, 1, 3, 1, -1);
    igraph_is_degree_sequence(&outseq, &inseq, &result);
    if (result) {
        return 5;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Directed degree sequence with different sums */
    igraph_vector_init_int_end(&outseq, -1, 0, 2, 3, 0, 4, 3, 1, 3, 4, 2, -1);
    igraph_vector_init_int_end(&inseq, -1, 0, 3, 1, 2, 2, 4, 4, 1, 3, 1, -1);
    igraph_is_degree_sequence(&outseq, &inseq, &result);
    if (result) {
        return 6;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /***** Testing igraph_is_graphical_degree_sequence *****/

    /* Valid undirected graphical degree sequence */
    igraph_vector_init_int_end(&outseq, -1, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    igraph_is_graphical_degree_sequence(&outseq, 0, &result);
    if (!result) {
        return 7;
    }
    igraph_vector_destroy(&outseq);

    /* Another valid undirected graphical degree sequence */
    igraph_vector_init_int_end(&outseq, -1, 4, 7, 4, 7, 7, 8, 9, 9, 4, 6, 5, -1);
    igraph_is_graphical_degree_sequence(&outseq, 0, &result);
    if (!result) {
        return 8;
    }
    igraph_vector_destroy(&outseq);

    /* Valid undirected degree sequence but not graphical */
    igraph_vector_init_int_end(&outseq, -1, 3, 3, -1);
    igraph_is_graphical_degree_sequence(&outseq, 0, &result);
    if (result) {
        return 9;
    }
    igraph_vector_destroy(&outseq);

    /* Valid directed graphical degree sequence */
    igraph_vector_init_int_end(&inseq, -1, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    igraph_vector_init_int_end(&outseq, -1, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (!result) {
        return 10;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Another valid directed graphical degree sequence */
    igraph_vector_init_int_end(&inseq, -1, 1, 3, 2, 1, 3, 4, 3, 3, 1, 3, -1);
    igraph_vector_init_int_end(&outseq, -1, 4, 1, 2, 3, 2, 3, 2, 3, 2, 2, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (!result) {
        return 11;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Yet another valid directed graphical degree sequence */
    igraph_vector_init_int_end(&inseq, -1, 7, 4, 6, 4, 7, 8, 8, 8, 7, 4, -1);
    igraph_vector_init_int_end(&outseq, -1, 8, 5, 6, 8, 6, 6, 5, 7, 5, 7, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (!result) {
        return 12;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Invalid directed graphical degree sequence when there is only one vertex
    * with a non-zero out-degree. Regression test for bug #851 */
    igraph_vector_init_int_end(&inseq, -1, 1, -1);
    igraph_vector_init_int_end(&outseq, -1, 1, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (result) {
        return 13;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Another invalid directed graphical degree sequence when there is only
     * one vertex with a non-zero out-degree. Regression test for bug #851 */
    igraph_vector_init_int_end(&inseq, -1, 2, 0, -1);
    igraph_vector_init_int_end(&outseq, -1, 0, 2, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (result) {
        return 14;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Another invalid directed graphical degree sequence when there is only
     * one vertex with a non-zero out-degree. Regression test for bug #851 */
    igraph_vector_init_int_end(&inseq, -1, 2, 2, -1);
    igraph_vector_init_int_end(&outseq, -1, 2, 2, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (result) {
        return 15;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    /* Valid directed graphical degree sequence. Regression test for bug #1092 */
    igraph_vector_init_int_end(&inseq, -1, 1, 0, 1, -1);
    igraph_vector_init_int_end(&outseq, -1, 0, 2, 0, -1);
    igraph_is_graphical_degree_sequence(&outseq, &inseq, &result);
    if (!result) {
        return 16;
    }
    igraph_is_graphical_degree_sequence(&inseq, &outseq, &result);
    if (!result) {
        return 17;
    }
    igraph_vector_destroy(&outseq);
    igraph_vector_destroy(&inseq);

    return 0;
}
