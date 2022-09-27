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
    igraph_vector_int_t v;
    igraph_t g;

    CHECK_ERROR(igraph_barabasi_game(&g, -10, /*power=*/ 1, 1, 0, 0, /*A=*/ 1, 0,
                               IGRAPH_BARABASI_BAG, /*start_from=*/ 0), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_barabasi_game(&g, 10, /*power=*/ 1, -2, 0, 0, /*A=*/ 1, 0,
                               IGRAPH_BARABASI_BAG, /*start_from=*/ 0), IGRAPH_EINVAL);
    igraph_vector_int_init(&v, 9);
    CHECK_ERROR(igraph_barabasi_game(&g, 10, /*power=*/ 1, 0, &v, 0, /*A=*/ 1, 0,
                               IGRAPH_BARABASI_BAG, /*start_from=*/ 0), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
