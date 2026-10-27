/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t g;
    igraph_matrix_t res, res2;
    igraph_vector_int_t roots, rootlevel;

    /* Two components: 0-1 and 2-3 */
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 2,3, -1);
    igraph_matrix_init(&res, 0, 0);
    igraph_matrix_init(&res2, 0, 0);

    /* 'roots' must not be modified, and repeated calls must give the same result. */
    igraph_vector_int_init_int(&roots, 2, 0, 2);
    igraph_vector_int_init_int(&rootlevel, 2, 1, 0);
    igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, &rootlevel);
    IGRAPH_ASSERT(VECTOR(roots)[0] == 0 && VECTOR(roots)[1] == 2);
    IGRAPH_ASSERT(MATRIX(res, 0, 1) == MATRIX(res, 2, 1) + 1);
    igraph_layout_reingold_tilford(&g, &res2, IGRAPH_ALL, &roots, &rootlevel);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res, &res2));
    igraph_vector_int_destroy(&roots);

    /* Invalid root IDs, including those of the vertices added internally. */
    igraph_vector_int_init_int(&roots, 1, 4);
    CHECK_ERROR(igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, NULL), IGRAPH_EINVVID);
    igraph_vector_int_destroy(&roots);

    igraph_vector_int_init_int(&roots, 2, 0, 4);
    CHECK_ERROR(igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, NULL), IGRAPH_EINVVID);
    CHECK_ERROR(igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, &rootlevel), IGRAPH_EINVVID);
    VECTOR(roots)[1] = -1;
    CHECK_ERROR(igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, NULL), IGRAPH_EINVVID);
    igraph_vector_int_destroy(&roots);

    /* Negative root levels. */
    igraph_vector_int_init_int(&roots, 2, 0, 2);
    VECTOR(rootlevel)[0] = 3; VECTOR(rootlevel)[1] = -1;
    CHECK_ERROR(igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, &rootlevel), IGRAPH_EINVAL);
    VECTOR(rootlevel)[0] = 2; VECTOR(rootlevel)[1] = -2;
    CHECK_ERROR(igraph_layout_reingold_tilford(&g, &res, IGRAPH_ALL, &roots, &rootlevel), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&roots);

    igraph_vector_int_destroy(&rootlevel);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
