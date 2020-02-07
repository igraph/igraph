/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

void print_matrix(igraph_matrix_t *m, FILE *f) {
    long int i, j;
    for (i = 0; i < igraph_matrix_nrow(m); i++) {
        for (j = 0; j < igraph_matrix_ncol(m); j++) {
            fprintf(f, " %.2f", MATRIX(*m, i, j));
        }
        fprintf(f, "\n");
    }
    fprintf(f, "==========\n");
}

int check_jaccard_all(const igraph_t* g, igraph_matrix_t* m,
                      igraph_neimode_t mode, igraph_bool_t loops) {
    igraph_vector_t pairs, res;
    long int i, j, k, n;
    igraph_eit_t eit;

    igraph_vector_init(&res, 0);

    /* First, query the similarities for all the vertices to a matrix */
    igraph_similarity_jaccard(g, m, igraph_vss_all(), mode, loops);

    /* Second, query the similarities for all pairs using a pair vector */
    n = igraph_vcount(g);
    igraph_vector_init(&pairs, 0);
    for (i = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--) {
            igraph_vector_push_back(&pairs, i);
            igraph_vector_push_back(&pairs, j);
        }
    }
    igraph_similarity_jaccard_pairs(g, &res, &pairs, mode, loops);
    for (i = 0, k = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--, k++) {
            if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
                fprintf(stderr, "Jaccard similarity calculation for vertex pair %ld-%ld "
                        "does not match the value in the full matrix (%.6f vs %.6f)\n",
                        i, j, VECTOR(res)[k], MATRIX(*m, i, j));
                return 1;
            }
        }
    }
    igraph_vector_destroy(&pairs);

    /* Third, query the similarities for all edges */
    igraph_similarity_jaccard_es(g, &res, igraph_ess_all(IGRAPH_EDGEORDER_FROM), mode, loops);
    igraph_eit_create(g, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &eit);
    k = 0;
    while (!IGRAPH_EIT_END(eit)) {
        long int eid = IGRAPH_EIT_GET(eit);
        i = IGRAPH_FROM(g, eid);
        j = IGRAPH_TO(g, eid);
        if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
            fprintf(stderr, "Jaccard similarity calculation for edge %ld-%ld (ID=%ld) "
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
    igraph_vector_t pairs, res;
    long int i, j, k, n;
    igraph_eit_t eit;

    igraph_vector_init(&res, 0);

    /* First, query the similarities for all the vertices to a matrix */
    igraph_similarity_dice(g, m, igraph_vss_all(), mode, loops);

    /* Second, query the similarities for all pairs using a pair vector */
    n = igraph_vcount(g);
    igraph_vector_init(&pairs, 0);
    for (i = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--) {
            igraph_vector_push_back(&pairs, i);
            igraph_vector_push_back(&pairs, j);
        }
    }
    igraph_similarity_dice_pairs(g, &res, &pairs, mode, loops);
    for (i = 0, k = 0; i < n; i++) {
        for (j = n - 1; j >= 0; j--, k++) {
            if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
                fprintf(stderr, "Dice similarity calculation for vertex pair %ld-%ld "
                        "does not match the value in the full matrix (%.6f vs %.6f)\n",
                        i, j, VECTOR(res)[k], MATRIX(*m, i, j));
                return 1;
            }
        }
    }
    igraph_vector_destroy(&pairs);

    /* Third, query the similarities for all edges */
    igraph_similarity_dice_es(g, &res, igraph_ess_all(IGRAPH_EDGEORDER_FROM), mode, loops);
    igraph_eit_create(g, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &eit);
    k = 0;
    while (!IGRAPH_EIT_END(eit)) {
        long int eid = IGRAPH_EIT_GET(eit);
        i = IGRAPH_FROM(g, eid);
        j = IGRAPH_TO(g, eid);
        if (fabs(VECTOR(res)[k] - MATRIX(*m, i, j)) > 1e-6) {
            fprintf(stderr, "Dice similarity calculation for edge %ld-%ld (ID=%ld) "
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

int main() {

    igraph_t g;
    igraph_matrix_t m;
    int ret;

    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 2, 1, 2, 0, 3, 0,
                 -1);

    igraph_matrix_init(&m, 0, 0);

    ret = check_jaccard_all(&g, &m, IGRAPH_ALL, 1);
    print_matrix(&m, stdout);
    if (ret) {
        return 1;
    }

    igraph_similarity_jaccard(&g, &m, igraph_vss_seq(1, 2), IGRAPH_ALL, 0);
    print_matrix(&m, stdout);

    ret = check_jaccard_all(&g, &m, IGRAPH_OUT, 1);
    print_matrix(&m, stdout);
    if (ret) {
        return 3;
    }

    ret = check_jaccard_all(&g, &m, IGRAPH_IN, 0);
    print_matrix(&m, stdout);
    if (ret) {
        return 4;
    }

    ret = check_dice_all(&g, &m, IGRAPH_ALL, 1);
    print_matrix(&m, stdout);
    if (ret) {
        return 5;
    }

    ret = check_dice_all(&g, &m, IGRAPH_OUT, 1);
    print_matrix(&m, stdout);
    if (ret) {
        return 6;
    }

    ret = check_dice_all(&g, &m, IGRAPH_IN, 0);
    print_matrix(&m, stdout);
    if (ret) {
        return 7;
    }

    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_ALL);
    print_matrix(&m, stdout);

    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_OUT);
    print_matrix(&m, stdout);

    igraph_similarity_inverse_log_weighted(&g, &m, igraph_vss_all(), IGRAPH_IN);
    print_matrix(&m, stdout);

    igraph_matrix_destroy(&m);
    igraph_destroy(&g);

    return 0;
}
