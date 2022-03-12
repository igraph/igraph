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
#include "test_utilities.inc"

int main() {
    igraph_vector_t alpha;
    igraph_matrix_t res;

    igraph_vector_init_int(&alpha, 2, 1, 1);
    igraph_matrix_init(&res, 0, 0);
    igraph_sample_dirichlet(0, &alpha, &res);
    igraph_matrix_print(&res);
    igraph_vector_destroy(&alpha);
    igraph_matrix_destroy(&res);

    VERIFY_FINALLY_STACK();
    igraph_vector_init(&alpha, 0);
    igraph_matrix_init(&res, 0, 0);

    CHECK_ERROR(igraph_sample_dirichlet(0, &alpha, &res), IGRAPH_EINVAL);
    igraph_vector_destroy(&alpha);
    igraph_vector_init_int(&alpha, 2, 1, 1);
    CHECK_ERROR(igraph_sample_dirichlet(-1, &alpha, &res), IGRAPH_EINVAL);

    igraph_vector_destroy(&alpha);
    igraph_matrix_destroy(&res);

    VERIFY_FINALLY_STACK();
    return 0;
}
