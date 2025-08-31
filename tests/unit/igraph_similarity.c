/* igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.h"

int check_jaccard_all(const igraph_t* g, igraph_matrix_t* m,
                      igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_vector_int_t pairs;
    igraph_vector_t res;
    igraph_int_t i, j, k, n;
    igraph_eit_t eit;

    igraph_vector_init(&res, 0);

    /* First, query the similarities for all the vertices to a matrix */
    igraph_similarity_jaccard(g, m, igraph_vss_all(), igraph_vss_all(), mode, loops);

    /* Second, query the similarities for all pairs using a pair vector */
    n = igraph_vcount(g);
    igraph_vector_int_init(&pairs, 0);
    for (i = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--) {
            igraph_vector_int_push_back(&pairs, i);
            igraph_vector_int_push_back(&pairs, j);
        }
    }
    igraph_similarity_jaccard_pairs(g, &res, &pairs, mode, loops);
    for (i = 0, k = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--, k++) {
            if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
                fprintf(stderr, "Jaccard similarity calculation for vertex pair %" IGRAPH_PRId "-%" IGRAPH_PRId " "
                        "does not match the value in the full matrix (%.6f vs %.6f)\n",
                        i, j, VECTOR(res)[k], MATRIX(*m, i, j));
                return 1;
            }
        }
    }
    igraph_vector_int_destroy(&pairs);

    /* Third, query the similarities for all edges */
    igraph_similarity_jaccard_es(g, &res, igraph_ess_all(IGRAPH_EDGEORDER_FROM), mode, loops);
    igraph_eit_create(g, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &eit);
    k = 0;
    while (!IGRAPH_EIT_END(eit)) {
        igraph_int_t eid = IGRAPH_EIT_GET(eit);
        i = IGRAPH_FROM(g, eid);
        j = IGRAPH_TO(g, eid);
        if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
            fprintf(stderr, "Jaccard similarity calculation for edge %" IGRAPH_PRId "-%" IGRAPH_PRId " (ID=%" IGRAPH_PRId ") "
                    "does not match the value in the full matrix (%.6f vs %.6f)\n",
                    i, j, eid, VECTOR(res)[k], MATRIX(*m, i, j));
            return 1;
        }
        IGRAPH_EIT_NEXT(eit);
        k++;
    }

    igraph_eit_destroy(&eit);

    igraph_vector_destroy(&res);

    return 0;
}

int check_dice_all(const igraph_t* g, igraph_matrix_t* m,
                   igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_vector_int_t pairs;
    igraph_vector_t res;
    igraph_int_t i, j, k, n;
    igraph_eit_t eit;

    igraph_vector_init(&res, 0);

    /* First, query the similarities for all the vertices to a matrix */
    igraph_similarity_dice(g, m, igraph_vss_all(), igraph_vss_all(), mode, loops);

    /* Second, query the similarities for all pairs using a pair vector */
    n = igraph_vcount(g);
    igraph_vector_int_init(&pairs, 0);
    for (i = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--) {
            igraph_vector_int_push_back(&pairs, i);
            igraph_vector_int_push_back(&pairs, j);
        }
    }
    igraph_similarity_dice_pairs(g, &res, &pairs, mode, loops);
    for (i = 0, k = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--, k++) {
            if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
                fprintf(stderr, "Dice similarity calculation for vertex pair %" IGRAPH_PRId "-%" IGRAPH_PRId " "
                        "does not match the value in the full matrix (%.6f vs %.6f)\n",
                        i, j, VECTOR(res)[k], MATRIX(*m, i, j));
                return 1;
            }
        }
    }
    igraph_vector_int_destroy(&pairs);

    /* Third, query the similarities for all edges */
    igraph_similarity_dice_es(g, &res, igraph_ess_all(IGRAPH_EDGEORDER_FROM), mode, loops);
    igraph_eit_create(g, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &eit);
    k = 0;
    while (!IGRAPH_EIT_END(eit)) {
        igraph_int_t eid = IGRAPH_EIT_GET(eit);
        i = IGRAPH_FROM(g, eid);
        j = IGRAPH_TO(g, eid);
        if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
            fprintf(stderr, "Dice similarity calculation for edge %" IGRAPH_PRId "-%" IGRAPH_PRId " (ID=%" IGRAPH_PRId ") "
                    "does not match the value in the full matrix (%.6f vs %.6f)\n",
                    i, j, eid, VECTOR(res)[k], MATRIX(*m, i, j));
            return 1;
        }
        IGRAPH_EIT_NEXT(eit);
        k++;
    }

    igraph_eit_destroy(&eit);

    igraph_vector_destroy(&res);

    return 0;
}

int main(void) {

    igraph_t g;
    igraph_matrix_t m;
    int ret;

    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 2, 1, 2, 0, 3, 0,
                 -1);

    igraph_matrix_init(&m, 0, 0);

    ret = check_jaccard_all(&g, &m, IGRAPH_ALL, 1);
    print_matrix(&m);
    if (ret) {
        return 1;
    }

    igraph_similarity_jaccard(&g, &m, igraph_vss_range(1, 3), igraph_vss_range(1, 3), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_matrix(&m);

    igraph_similarity_jaccard(&g, &m, igraph_vss_range(1, 3), igraph_vss_1(1), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_matrix(&m);

    ret = check_jaccard_all(&g, &m, IGRAPH_OUT, 1);
    print_matrix(&m);
    if (ret) {
        return 3;
    }

    ret = check_jaccard_all(&g, &m, IGRAPH_IN, 0);
    print_matrix(&m);
    if (ret) {
        return 4;
    }

    ret = check_dice_all(&g, &m, IGRAPH_ALL, 1);
    print_matrix(&m);
    if (ret) {
        return 5;
    }

    ret = check_dice_all(&g, &m, IGRAPH_OUT, 1);
    print_matrix(&m);
    if (ret) {
        return 6;
    }

    ret = check_dice_all(&g, &m, IGRAPH_IN, 0);
    print_matrix(&m);
    if (ret) {
        return 7;
    }

    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_ALL);
    print_matrix(&m);

    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_OUT);
    print_matrix(&m);

    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_IN);
    print_matrix(&m);

    igraph_matrix_destroy(&m);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
