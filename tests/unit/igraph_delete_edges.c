/* igraph library.
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
    igraph_es_t es;

    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, -1);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 3, 2, -1);

    /* error test, no such edge to delete */
    CHECK_ERROR(igraph_delete_edges(&g, es), IGRAPH_EINVAL);
    if (igraph_ecount(&g) != 3) {
        return 3;
    }

    /* error test, invalid vertex ID */
    igraph_es_destroy(&es);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 10, 2, -1);
    CHECK_ERROR(igraph_delete_edges(&g, es), IGRAPH_EINVVID);
    if (igraph_ecount(&g) != 3) {
        return 5;
    }

    /* error test, invalid (odd) length */
    igraph_es_destroy(&es);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 0, 1, 2, -1);
    CHECK_ERROR(igraph_delete_edges(&g, es), IGRAPH_EINVAL);
    if (igraph_ecount(&g) != 3) {
        return 7;
    }

    igraph_es_destroy(&es);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
