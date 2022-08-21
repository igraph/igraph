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
    igraph_t pattern, target;
    igraph_bool_t iso;
    igraph_vector_int_t map;
    igraph_vector_int_list_t maps;

    igraph_small(&target, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 6,
                 1, 4, 1, 2,
                 2, 3,
                 3, 4, 3, 5, 3, 7, 3, 8,
                 4, 5, 4, 6,
                 5, 6, 5, 8,
                 7, 8,
                 -1);
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_list_init(&maps, 0);
    igraph_small(&pattern, 0, IGRAPH_DIRECTED, -1);
    /* pattern and target differ in directedness */
    CHECK_ERROR(igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ 0,
                                      &iso, &map, &maps, /*induced=*/ 0,
                                      /*time_limit=*/ 0), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_list_destroy(&maps);
    igraph_destroy(&pattern);
    igraph_destroy(&target);
    VERIFY_FINALLY_STACK();

    return 0;
}
