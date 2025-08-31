/*
   igraph library.
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

int main(void) {
    igraph_t g;
    igraph_sparsemat_t res;

    igraph_sparsemat_init(&res, 0, 0, 0);

    printf("Graph with no vertices:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
    igraph_get_stochastic_sparse(&g, &res, /* column_wise = */ 0, /* weights = */ NULL);
    igraph_sparsemat_print(&res, stdout);
    igraph_destroy(&g);

    printf("\nGraph with 4 vertices, rowwise:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0,1, 1,3, 1,3, 2,2, 2,2, 2,3, 3,0, -1);
    igraph_get_stochastic_sparse(&g, &res, /* column_wise = */ 0, /* weights = */ NULL);
    igraph_sparsemat_print(&res, stdout);
    igraph_destroy(&g);

    printf("\nColumnwise:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0,1, 1,3, 1,3, 2,2, 2,2, 2,3, 3,0, -1);
    igraph_get_stochastic_sparse(&g, &res, /* column_wise = */ 1, /* weights = */ NULL);
    igraph_sparsemat_print(&res, stdout);
    igraph_destroy(&g);

    igraph_sparsemat_destroy(&res);

    VERIFY_FINALLY_STACK();
    return 0;
}
