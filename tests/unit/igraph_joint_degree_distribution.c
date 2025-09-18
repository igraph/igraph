/*
   igraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

/* Check vector equality with tolerances. Consider NaN values equal. */
igraph_bool_t vector_eq(const igraph_vector_t *a, const igraph_vector_t *b) {
    igraph_int_t na = igraph_vector_size(a);
    igraph_int_t nb = igraph_vector_size(b);
    if (na != nb) {
        return false;
    }
    for (igraph_int_t i=0; i < na; i++) {
        if (isnan(VECTOR(*a)[i]) && isnan(VECTOR(*b)[i])) {
            continue;
        }
        if (! igraph_almost_equals(VECTOR(*a)[i], VECTOR(*b)[i], 1e-12)) {
            return false;
        }
    }
    return true;
}

/* Compare to igraph_joint_degree_matrix() */
void check_jdm(const igraph_t *g, const igraph_vector_t *weights) {
    igraph_matrix_t jdm, p;
    igraph_int_t vcount = igraph_vcount(g);
    igraph_int_t nrow, ncol, n;

    igraph_matrix_init(&jdm, 0, 0);
    igraph_matrix_init(&p, 0, 0);

    igraph_joint_degree_matrix(g, weights, &jdm, -1, -1);
    igraph_joint_degree_distribution(g, weights, &p, IGRAPH_OUT, IGRAPH_IN, true, /*normalized*/ false, -1, -1);

    nrow = igraph_matrix_nrow(&p);
    ncol = igraph_matrix_ncol(&p);

    if (vcount > 0) {
        igraph_vector_t v;

        IGRAPH_ASSERT(nrow > 0);
        IGRAPH_ASSERT(ncol > 0);
        nrow--;
        ncol--;

        IGRAPH_ASSERT(igraph_matrix_nrow(&jdm) == nrow);
        IGRAPH_ASSERT(igraph_matrix_ncol(&jdm) == ncol);

        igraph_vector_init(&v, 0);

        igraph_matrix_get_col(&p, &v, 0);
        IGRAPH_ASSERT(igraph_vector_isnull(&v));

        igraph_matrix_get_row(&p, &v, 0);
        IGRAPH_ASSERT(igraph_vector_isnull(&v));

        igraph_vector_destroy(&v);

        igraph_matrix_remove_row(&p, 0);
        igraph_matrix_remove_col(&p, 0);
    } else {
        // vcount == 0
        IGRAPH_ASSERT(nrow == 0);
        IGRAPH_ASSERT(ncol == 0);
    }

    n = nrow < ncol ? nrow : ncol;

    if (! igraph_is_directed(g)) {
        for (igraph_int_t i=0; i < n; i++) {
            MATRIX(jdm, i, i) *= 2;
        }
    }

    igraph_real_t total;
    if (weights) {
        total = igraph_vector_sum(weights);
    } else {
        total = igraph_ecount(g);
    }
    if (! igraph_is_directed(g)) {
        total *= 2;
    }

    IGRAPH_ASSERT(igraph_matrix_sum(&p) == total);

    IGRAPH_ASSERT(igraph_matrix_all_e(&jdm, &p));

    igraph_matrix_destroy(&p);
    igraph_matrix_destroy(&jdm);
}

void check_assort_i(const igraph_t *g, const igraph_vector_t *weights,
                    igraph_neimode_t from_mode, igraph_neimode_t to_mode) {
    igraph_matrix_t p;
    igraph_real_t r1, r2;
    igraph_int_t nrow, ncol;
    igraph_vector_t dfrom, dto;
    igraph_vector_t a, b;
    igraph_bool_t directed = igraph_is_directed(g);

    igraph_vector_init(&dfrom, 0);
    igraph_vector_init(&dto, 0);
    igraph_vector_init(&a, 0);
    igraph_vector_init(&b, 0);

    igraph_strength(g, &dfrom, igraph_vss_all(), from_mode, IGRAPH_LOOPS, NULL);
    igraph_strength(g, &dto, igraph_vss_all(), to_mode, IGRAPH_LOOPS, NULL);

    igraph_matrix_init(&p, 0, 0);
    igraph_joint_degree_distribution(g, weights, &p, from_mode, to_mode, true, /*normalized*/ true, -1, -1);

    nrow = igraph_matrix_nrow(&p);
    ncol = igraph_matrix_ncol(&p);

    igraph_assortativity(g, weights, &dfrom, directed ? &dto : NULL, &r1, IGRAPH_DIRECTED, /*normalized*/ false);

    igraph_matrix_rowsum(&p, &a);
    igraph_matrix_colsum(&p, &b);

    r2 = 0;
    for (igraph_int_t i=0; i < nrow; i++) {
        for (igraph_int_t j=0; j < ncol; j++) {
            r2 += (MATRIX(p, i, j) - VECTOR(a)[i] * VECTOR(b)[j]) * (igraph_real_t) i * (igraph_real_t) j;
        }
    }
    // printf("Assortativity: %g == %g\n", r1, r2);
    IGRAPH_ASSERT(igraph_almost_equals(r1, r2, 1e-13));

    igraph_matrix_destroy(&p);

    igraph_vector_destroy(&b);
    igraph_vector_destroy(&a);
    igraph_vector_destroy(&dto);
    igraph_vector_destroy(&dfrom);
}

/* Compare to igraph_assortativity() */
void check_assort(const igraph_t *g, const igraph_vector_t *weights) {
    if (igraph_is_directed(g)) {
        check_assort_i(g, weights, IGRAPH_OUT, IGRAPH_IN);
        check_assort_i(g, weights, IGRAPH_IN, IGRAPH_OUT);
        check_assort_i(g, weights, IGRAPH_OUT, IGRAPH_OUT);
        check_assort_i(g, weights, IGRAPH_IN, IGRAPH_IN);
        check_assort_i(g, weights, IGRAPH_ALL, IGRAPH_IN);
        check_assort_i(g, weights, IGRAPH_ALL, IGRAPH_OUT);
        check_assort_i(g, weights, IGRAPH_OUT, IGRAPH_ALL);
        check_assort_i(g, weights, IGRAPH_IN, IGRAPH_ALL);
    } else {
        check_assort_i(g, weights, IGRAPH_ALL, IGRAPH_ALL);
    }
}

void check_knnk_i(const igraph_t *g, const igraph_vector_t *weights, igraph_neimode_t from_mode, igraph_neimode_t to_mode) {
    igraph_matrix_t p;
    igraph_int_t nrow, ncol;
    igraph_vector_t knnk, knnk2;
    igraph_vector_t q;

    igraph_vector_init(&knnk, 0);
    igraph_vector_init(&q, 0);

    igraph_matrix_init(&p, 0, 0);
    igraph_joint_degree_distribution(g, weights, &p, from_mode, to_mode, true, /*normalized*/ true, -1, -1);

    nrow = igraph_matrix_nrow(&p);
    ncol = igraph_matrix_ncol(&p);

    igraph_degree_correlation_vector(g, weights, &knnk, from_mode, to_mode, /*directed_neighbors*/ true);
    IGRAPH_ASSERT(igraph_vector_size(&knnk) == nrow);

    igraph_vector_init(&knnk2, nrow);
    igraph_matrix_rowsum(&p, &q);
    for (igraph_int_t k=0; k < nrow; k++) {
        for (igraph_int_t j=0; j < ncol; j++) {
            VECTOR(knnk2)[k] += j * MATRIX(p, k, j);
        }
    }

    igraph_vector_div(&knnk2, &q);

    /*
    printf("%d - %d\n", from_mode, to_mode);
    print_vector(&knnk);
    print_vector(&knnk2);
    */
    IGRAPH_ASSERT(vector_eq(&knnk, &knnk2));

    igraph_vector_destroy(&knnk2);

    igraph_matrix_destroy(&p);

    igraph_vector_destroy(&q);
    igraph_vector_destroy(&knnk);
}

/* Compare to igraph_degree_correlation_vector() */
void check_knnk(const igraph_t *g, const igraph_vector_t *weights) {
    if (igraph_is_directed(g)) {
        check_knnk_i(g, weights, IGRAPH_OUT, IGRAPH_IN);
        check_knnk_i(g, weights, IGRAPH_IN, IGRAPH_OUT);
        check_knnk_i(g, weights, IGRAPH_OUT, IGRAPH_OUT);
        check_knnk_i(g, weights, IGRAPH_IN, IGRAPH_IN);
        check_knnk_i(g, weights, IGRAPH_ALL, IGRAPH_IN);
        check_knnk_i(g, weights, IGRAPH_ALL, IGRAPH_OUT);
        check_knnk_i(g, weights, IGRAPH_OUT, IGRAPH_ALL);
        check_knnk_i(g, weights, IGRAPH_IN, IGRAPH_ALL);
    } else {
        check_knnk_i(g, weights, IGRAPH_ALL, IGRAPH_ALL);
    }
}

int main(void) {
    igraph_t dg, ug;
    igraph_vector_t weights;
    igraph_matrix_t p;

    igraph_rng_seed(igraph_rng_default(), 137);

    igraph_matrix_init(&p, 0, 0);

    /* Directed */

    igraph_empty(&dg, 0, IGRAPH_DIRECTED);
    check_jdm(&dg, NULL);
    igraph_destroy(&dg);

    igraph_empty(&dg, 1, IGRAPH_DIRECTED);
    check_jdm(&dg, NULL);
    igraph_destroy(&dg);

    igraph_small(&dg, 2, IGRAPH_DIRECTED, 0,0, 0,1, 0,1, 1,0, -1);
    check_jdm(&dg, NULL);
    check_assort(&dg, NULL);
    check_knnk(&dg, NULL);
    igraph_destroy(&dg);

    igraph_erdos_renyi_game_gnm(&dg, 10, 30, IGRAPH_DIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    check_jdm(&dg, NULL);
    check_assort(&dg, NULL);
    check_knnk(&dg, NULL);

    igraph_vector_init_range(&weights, 0, igraph_ecount(&dg));
    check_jdm(&dg, &weights);
    check_assort(&dg, &weights);
    check_knnk(&dg, &weights);
    igraph_vector_destroy(&weights);

    igraph_joint_degree_distribution(&dg, NULL, &p, IGRAPH_ALL, IGRAPH_ALL, false, false, -1, -1);
    print_matrix(&p);

    igraph_joint_degree_distribution(&dg, NULL, &p, IGRAPH_ALL, IGRAPH_ALL, false, false, 3, 10);
    print_matrix(&p);

    igraph_destroy(&dg);

    /* Undirected */

    igraph_empty(&ug, 0, IGRAPH_UNDIRECTED);
    check_jdm(&ug, NULL);
    igraph_destroy(&ug);

    igraph_empty(&ug, 1, IGRAPH_UNDIRECTED);
    check_jdm(&ug, NULL);
    igraph_destroy(&ug);

    igraph_small(&ug, 2, IGRAPH_UNDIRECTED, 0,1, 0,1, 1,1, -1);
    check_jdm(&ug, NULL);
    check_assort(&ug, NULL);
    check_knnk(&ug, NULL);
    igraph_destroy(&ug);

    igraph_erdos_renyi_game_gnm(&ug, 10, 30, IGRAPH_UNDIRECTED, IGRAPH_LOOPS_SW, IGRAPH_EDGE_UNLABELED);
    check_jdm(&ug, NULL);
    check_assort(&ug, NULL);
    check_knnk(&ug, NULL);

    igraph_vector_init_range(&weights, 0, igraph_ecount(&ug));
    check_jdm(&ug, &weights);
    check_assort(&ug, &weights);
    check_knnk(&ug, &weights);
    igraph_vector_destroy(&weights);

    igraph_joint_degree_distribution(&ug, NULL, &p, IGRAPH_ALL, IGRAPH_ALL, false, false, -1, -1);
    print_matrix(&p);

    igraph_joint_degree_distribution(&ug, NULL, &p, IGRAPH_ALL, IGRAPH_ALL, false, false, 3, 10);
    print_matrix(&p);

    igraph_destroy(&ug);

    igraph_matrix_destroy(&p);

    VERIFY_FINALLY_STACK();

    return 0;
}
