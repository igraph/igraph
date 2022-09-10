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

int main() {
    igraph_t g;
    igraph_vector_bool_t types;
    igraph_vector_t weights;

    igraph_integer_t matching_size;
    igraph_real_t weighted_size;

    igraph_vector_int_t matching;

    igraph_small(&g, 4, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3,
                 -1);

    igraph_vector_bool_init(&types, 4);
    VECTOR(types)[0] = 0;
    VECTOR(types)[1] = 1;
    VECTOR(types)[2] = 0;
    VECTOR(types)[3] = 1;

    igraph_vector_int_init(&matching, 0);

    igraph_vector_init(&weights, igraph_vcount(&g));
    igraph_vector_fill(&weights, 1.0);

    // Test incorrect types
    CHECK_ERROR(igraph_maximum_bipartite_matching(&g, &types, &matching_size, NULL, &matching, NULL, 0), IGRAPH_EINVAL);

    // Test incorrect types for weighted graph
    VECTOR(types)[2] = 0;
    CHECK_ERROR(igraph_maximum_bipartite_matching(&g, &types, &matching_size, &weighted_size, &matching, &weights, 0), IGRAPH_EINVAL);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&matching);

    igraph_vector_bool_destroy(&types);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
