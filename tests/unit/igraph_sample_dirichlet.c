/* IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void check_result(igraph_matrix_t *res, igraph_integer_t n)
{
    igraph_vector_t colsum;
    igraph_vector_init(&colsum, 0);

    IGRAPH_ASSERT(igraph_matrix_min(res) >= 0);
    IGRAPH_ASSERT(igraph_matrix_max(res) <= 1.0);
    igraph_matrix_colsum(res, &colsum);
    IGRAPH_ASSERT(igraph_vector_size(&colsum) == n);
    for (igraph_integer_t i = 0; i < igraph_vector_size(&colsum); i++) {
        IGRAPH_ASSERT(igraph_almost_equals(VECTOR(colsum)[i], 1, 0.0000001));
    }
    igraph_vector_destroy(&colsum);
}


int main(void) {
    igraph_vector_t alpha;
    igraph_matrix_t res;

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_matrix_init(&res, 0, 0);

    printf("Zero vectors to sample should return empty matrix:\n");
    igraph_vector_init_int(&alpha, 2, 1, 1);
    igraph_sample_dirichlet(0, &alpha, &res);
    igraph_matrix_print(&res);
    igraph_vector_destroy(&alpha);

    printf("Check if result vectors add up to one.\n");
    igraph_vector_init_int(&alpha, 5, 1, 2, 3, 4, 5);
    igraph_sample_dirichlet(100, &alpha, &res);
    check_result(&res, 100);
    igraph_vector_destroy(&alpha);

    printf("Distribution localized at 0.5, 0.5:\n");
    igraph_vector_init_real(&alpha, 2, 1e30, 1e30);
    igraph_sample_dirichlet(2, &alpha, &res);
    igraph_matrix_print(&res);
    igraph_vector_destroy(&alpha);

    VERIFY_FINALLY_STACK();

    printf("Check if too short parameter vector is handled correctly.\n");
    igraph_vector_init(&alpha, 0);
    CHECK_ERROR(igraph_sample_dirichlet(0, &alpha, &res), IGRAPH_EINVAL);
    igraph_vector_destroy(&alpha);

    printf("Check if negative number of samples is handled correctly.\n");
    igraph_vector_init_int(&alpha, 2, 1, 1);
    CHECK_ERROR(igraph_sample_dirichlet(-1, &alpha, &res), IGRAPH_EINVAL);
    igraph_vector_destroy(&alpha);

    printf("Check if negative alpha is handled correctly.\n");
    igraph_vector_init_int(&alpha, 2, -1, 1);
    CHECK_ERROR(igraph_sample_dirichlet(0, &alpha, &res), IGRAPH_EINVAL);
    igraph_vector_destroy(&alpha);

    igraph_matrix_destroy(&res);

    VERIFY_FINALLY_STACK();
    return 0;
}
