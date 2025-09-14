/*
   igraph library.
   Copyright (C) 2022 The igraph development team

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

#include "test_utilities.h"

void print_cliques_and_thresholds(igraph_vector_int_list_t* cliques, igraph_vector_t* thresholds) {
    printf("Cliques:\n");
    print_vector_int_list(cliques);
    printf("Thresholds:\n");
    print_vector(thresholds);
}

void print_cliques_and_mu(igraph_vector_int_list_t* cliques, igraph_vector_t* mu) {
    printf("Cliques:\n");
    print_vector_int_list(cliques);
    printf("Mu:\n");
    print_vector_format(mu, stdout, "%.5f");
}

void test_graphlets_candidate_basis_simple(void) {
    igraph_t g;
    igraph_vector_int_list_t cliques;
    igraph_vector_t thresholds;
    igraph_vector_t weights;
    igraph_int_t eid;

    igraph_full(&g, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_int_list_init(&cliques, 0);
    igraph_vector_init(&thresholds, 0);
    igraph_vector_init(&weights, igraph_ecount(&g));
    igraph_vector_fill(&weights, 1);

    igraph_graphlets_candidate_basis(&g, &weights, &cliques, &thresholds);
    print_cliques_and_thresholds(&cliques, &thresholds);

    igraph_get_eid(&g, &eid, 0, 1, IGRAPH_UNDIRECTED, /* error = */ 0);
    VECTOR(weights)[eid] = 2;

    igraph_graphlets_candidate_basis(&g, &weights, &cliques, &thresholds);
    print_cliques_and_thresholds(&cliques, &thresholds);

    igraph_vector_destroy(&weights);
    igraph_vector_int_list_destroy(&cliques);
    igraph_vector_destroy(&thresholds);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
}

void test_graphlets_filtering(void) {
    igraph_t g;
    igraph_vector_int_list_t cliques;
    igraph_vector_t thresholds;
    igraph_vector_t weights_vec;
    igraph_real_t weights[] = { 8, 8, 8, 5, 5, 5, 5, 5 };

    igraph_small(&g, 5, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4, -1);
    igraph_vector_int_list_init(&cliques, 0);
    igraph_vector_init(&thresholds, 0);
    weights_vec = igraph_vector_view(weights, sizeof(weights) / sizeof(weights[0]));

    igraph_graphlets_candidate_basis(&g, &weights_vec, &cliques, &thresholds);
    print_cliques_and_thresholds(&cliques, &thresholds);

    igraph_vector_int_list_destroy(&cliques);
    igraph_vector_destroy(&thresholds);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
}

void test_zachary_random_weights(void) {
    igraph_t g;
    igraph_vector_int_list_t cliques;
    igraph_vector_t thresholds;
    igraph_vector_t weights_vec;
    igraph_real_t weights[] = {
        1, 5, 1, 1, 2, 4, 2, 2, 1, 4, 1, 5, 4, 2, 2, 3, 1, 1, 3, 4, 5, 5, 5, 4,
        2, 4, 3, 2, 1, 2, 3, 2, 4, 4, 2, 5, 4, 5, 4, 2, 2, 3, 1, 5, 2, 2, 2, 4,
        3, 5, 2, 2, 2, 5, 1, 1, 4, 5, 2, 1, 5, 4, 4, 1, 3, 3, 5, 5, 4, 5, 4, 2,
        2, 1, 2, 5, 5, 4
    };

    igraph_famous(&g, "zachary");
    igraph_vector_int_list_init(&cliques, 0);
    igraph_vector_init(&thresholds, 0);
    weights_vec = igraph_vector_view(weights, sizeof(weights) / sizeof(weights[0]));

    igraph_graphlets_candidate_basis(&g, &weights_vec, &cliques, &thresholds);
    print_cliques_and_thresholds(&cliques, &thresholds);

    igraph_vector_int_list_destroy(&cliques);
    igraph_vector_destroy(&thresholds);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

}

void test_projection(void) {
    igraph_t g;
    igraph_vector_int_list_t cliques;
    igraph_vector_t mu;
    igraph_vector_t weights_vec;
    igraph_real_t weights[] = { 2, 2, 3, 1, 1, 4, 4, 4 };

    igraph_small(&g, 5, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4, -1);
    igraph_vector_int_list_init(&cliques, 0);
    igraph_vector_init(&mu, 0);
    weights_vec = igraph_vector_view(weights, sizeof(weights) / sizeof(weights[0]));

    igraph_graphlets(&g, &weights_vec, &cliques, &mu, 1000);
    print_cliques_and_mu(&cliques, &mu);

    igraph_vector_int_list_destroy(&cliques);
    igraph_vector_destroy(&mu);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
}

int main(void) {
    test_graphlets_candidate_basis_simple();
    test_graphlets_filtering();
    test_zachary_random_weights();
    test_projection();

    return 0;
}
