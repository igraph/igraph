/* IGraph library.
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
    igraph_integer_t edges3[] = {0, 1, 1, 2, 3, 4, 5, 6, 6, 5, 2, 4, 1, 6, 0, 3 };
    igraph_vector_int_t edges;
    igraph_vector_bool_t types;
    igraph_t g;
    igraph_integer_t i;

    igraph_vector_int_view(&edges, edges3, sizeof(edges3) / sizeof(edges3[0]));
    igraph_vector_bool_init(&types, igraph_vector_int_max(&edges) + 1);
    for (i = 0; i < igraph_vector_bool_size(&types); i++) {
        VECTOR(types)[i] = i % 2;
    }
    /* Not a bipartite graph, vertices 2 and 4 are connected */
    CHECK_ERROR(igraph_create_bipartite(&g, &types, &edges, /*directed=*/ 1), IGRAPH_EINVAL);

    igraph_vector_bool_destroy(&types);
    VERIFY_FINALLY_STACK();

    return 0;
}
