/*
   igraph library.
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

int main(void) {

    igraph_t g;
    igraph_vector_int_t v;
    igraph_vector_int_t v2, v3;

    /* Initialize the library. */
    igraph_setup();

    igraph_barabasi_game(&g, 10, /*power=*/ 1, 2, 0, 0, /*A=*/ 1, 1,
                         IGRAPH_BARABASI_BAG, /*start_from=*/ 0);
    if (igraph_ecount(&g) != 18) {
        return 1;
    }
    if (igraph_vcount(&g) != 10) {
        return 2;
    }
    if (!igraph_is_directed(&g)) {
        return 3;
    }

    igraph_vector_int_init(&v, 0);
    igraph_get_edgelist(&g, &v, 0);
    for (igraph_int_t i = 0; i < igraph_ecount(&g); i++) {
        if (VECTOR(v)[2 * i] <= VECTOR(v)[2 * i + 1]) {
            return 4;
        }
    }
    igraph_vector_int_destroy(&v);
    igraph_destroy(&g);

    /* out-degree sequence */
    igraph_vector_int_init_int(&v3, 10, 0, 1, 3, 3, 4, 5, 6, 7, 8, 9);

    igraph_barabasi_game(&g, 10, /*power=*/ 1, 0, &v3, 0, /*A=*/ 1, 1,
                         IGRAPH_BARABASI_BAG, /*start_from=*/ 0);
    if (igraph_ecount(&g) != igraph_vector_int_sum(&v3)) {
        return 5;
    }
    igraph_vector_int_init(&v2, 0);
    igraph_degree(&g, &v2, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    for (igraph_int_t i = 0; i < igraph_vcount(&g); i++) {
        if (VECTOR(v3)[i] != VECTOR(v2)[i]) {
            igraph_vector_int_print(&v3);
            printf("\n");
            igraph_vector_int_print(&v2);
            return 6;
        }
    }
    igraph_vector_int_destroy(&v3);
    igraph_vector_int_destroy(&v2);
    igraph_destroy(&g);

    /* outpref, we cannot really test this quantitatively,
       would need to set random seed */
    igraph_barabasi_game(&g, 10, /*power=*/ 1, 2, 0, 1, /*A=*/ 1, 1,
                         IGRAPH_BARABASI_BAG, /*start_from=*/ 0);
    igraph_vector_int_init(&v, 0);
    igraph_get_edgelist(&g, &v, 0);
    for (igraph_int_t i = 0; i < igraph_ecount(&g); i++) {
        if (VECTOR(v)[2 * i] <= VECTOR(v)[2 * i + 1]) {
            return 7;
        }
    }
    if (!igraph_is_directed(&g)) {
        return 8;
    }
    igraph_vector_int_destroy(&v);
    igraph_destroy(&g);

    return 0;
}
