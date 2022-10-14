/*
   IGraph library.
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

void call_and_print(igraph_integer_t n, igraph_vector_int_t *mlist, igraph_vector_list_t *rholist, igraph_matrix_list_t *pref_matrix_list, igraph_real_t p) {
    igraph_t result;
    igraph_hsbm_list_game(&result, n, mlist, rholist, pref_matrix_list, p);
    print_graph_canon(&result);
    printf("\n");
    igraph_destroy(&result);
}


int main(void) {
    igraph_matrix_t pref_matrix;
    igraph_matrix_list_t pref_matrix_list;
    igraph_vector_t rho;
    igraph_vector_list_t rholist;
    igraph_vector_int_t mlist;

    igraph_vector_list_init(&rholist, 0);
    igraph_vector_init_int(&rho, 1, 1);
    igraph_vector_list_push_back(&rholist, &rho);

    igraph_matrix_list_init(&pref_matrix_list, 0);
    igraph_matrix_init(&pref_matrix, 1, 1);
    MATRIX(pref_matrix, 0, 0) = 1;
    igraph_matrix_list_push_back(&pref_matrix_list, &pref_matrix);

    igraph_vector_int_init_int(&mlist, 1, 1);
    printf("One block, one vertex.\n");
    call_and_print(1, &mlist, &rholist, &pref_matrix_list, 0);

    igraph_vector_list_clear(&rholist);
    igraph_matrix_list_clear(&pref_matrix_list);
    igraph_vector_int_destroy(&mlist);

    {
        igraph_vector_init_real(&rho, 3, 0.6, 0.4, 0.0);
        igraph_vector_list_push_back(&rholist, &rho);

        igraph_real_t elems[] = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        matrix_init_real_row_major(&pref_matrix, 3, 3, elems);
        igraph_matrix_list_push_back(&pref_matrix_list, &pref_matrix);

        igraph_vector_int_init_int(&mlist, 1, 10);
        printf("One block, two clusters, 6 and 4 vertices in cluster, complete bipartite.\n");
        call_and_print(10, &mlist, &rholist, &pref_matrix_list, 0);

        igraph_vector_list_clear(&rholist);
        igraph_matrix_list_clear(&pref_matrix_list);
        igraph_vector_int_destroy(&mlist);
    }

    {
        igraph_vector_init_real(&rho, 3, 0.6, 0.4, 0.0);
        igraph_vector_list_push_back(&rholist, &rho);
        igraph_vector_init_real(&rho, 3, 0.6, 0.4, 0.0);
        igraph_vector_list_push_back(&rholist, &rho);

        igraph_real_t elems[] = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        matrix_init_real_row_major(&pref_matrix, 3, 3, elems);
        igraph_matrix_list_push_back(&pref_matrix_list, &pref_matrix);
        matrix_init_real_row_major(&pref_matrix, 3, 3, elems);
        igraph_matrix_list_push_back(&pref_matrix_list, &pref_matrix);

        igraph_vector_int_init_int(&mlist, 2, 5, 5);
        printf("Two blocks, two clusters each, 3 and 2 vertices in cluster, each vertex connected to every other, except those in the same cluster.\n");
        call_and_print(10, &mlist, &rholist, &pref_matrix_list, 1);

        igraph_vector_list_destroy(&rholist);
        igraph_matrix_list_destroy(&pref_matrix_list);
        igraph_vector_int_destroy(&mlist);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
