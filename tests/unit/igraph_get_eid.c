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
    igraph_t g;
    igraph_integer_t eid;

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);

    /* NON-EXISTENT EDGE */
    CHECK_ERROR(igraph_get_eid(&g, &eid, 5, 6, IGRAPH_UNDIRECTED, /*error=*/ 1), IGRAPH_EINVAL);

    /* INVALID VERTEX ID */
    CHECK_ERROR(igraph_get_eid(&g, &eid, 171, 6, IGRAPH_UNDIRECTED, /*error=*/ 1), IGRAPH_EINVVID);

    /* INVALID VERTEX ID even if error == 0 */
    CHECK_ERROR(igraph_get_eid(&g, &eid, 171, 6, IGRAPH_UNDIRECTED, /*error=*/ 0), IGRAPH_EINVVID);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
