/*
   igraph library.
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

int main(void) {

    igraph_t g;
    igraph_matrix_t m;
    igraph_vector_int_t pairs;
    igraph_vector_t res;
    igraph_int_t i, j, n;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 2, 1, 2, 0, 3, 0,
                 -1);

    igraph_matrix_init(&m, 0, 0);
    igraph_vector_init(&res, 0);
    igraph_vector_int_init(&pairs, 0);

    n = igraph_vcount(&g);
    for (i = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--) {
            igraph_vector_int_push_back(&pairs, i);
            igraph_vector_int_push_back(&pairs, j);
        }
    }

    printf("Jaccard similarity:\n");
    igraph_similarity_jaccard(&g, &m, igraph_vss_range(1, 3), igraph_vss_range(1, 3), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_matrix_printf(&m, "%.2f");

    printf("\nJaccard similarity, pairs:\n");
    igraph_similarity_jaccard_pairs(&g, &res, &pairs, IGRAPH_ALL, 0);
    igraph_vector_print(&res);

    printf("\nJaccard similarity with edge selector:\n");
    igraph_similarity_jaccard_es(&g, &res, igraph_ess_all(IGRAPH_EDGEORDER_FROM), IGRAPH_IN, 0);
    igraph_vector_print(&res);

    printf("\nDice similarity:\n");
    igraph_similarity_dice(&g, &m, igraph_vss_range(1, 3), igraph_vss_range(1, 3), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_matrix_printf(&m, "%.2f");

    printf("\nDice similarity, pairs:\n");
    igraph_similarity_dice_pairs(&g, &res, &pairs, IGRAPH_ALL, 0);
    igraph_vector_print(&res);

    printf("\nDice similarity with edge selector:\n");
    igraph_similarity_dice_es(&g, &res, igraph_ess_all(IGRAPH_EDGEORDER_FROM), IGRAPH_IN, 0);
    igraph_vector_print(&res);

    printf("\nWeighted inverse log similarity:\n");
    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_ALL);
    igraph_matrix_printf(&m, "%.2f");

    igraph_matrix_destroy(&m);
    igraph_destroy(&g);
    igraph_vector_destroy(&res);
    igraph_vector_int_destroy(&pairs);

    return 0;
}
