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

    igraph_vector_int_init(&v, 9);
    /* error: IGRAPH_EINVEVECTOR */
    VECTOR(v)[8] = 0;
    CHECK_ERROR(igraph_create(&g, &v, 0, 0), IGRAPH_EINVEVECTOR);

    /* error: IGRAPH_EINVVID */
    igraph_vector_int_resize(&v, 8);
    VECTOR(v)[7] = -1;
    CHECK_ERROR(igraph_create(&g, &v, 10, 1), IGRAPH_EINVVID);

    igraph_vector_int_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
