/*
   IGraph library.
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
    igraph_t g1, g2;
    igraph_bool_t res;
    int err;

    /* undirected graphs */
    igraph_small(&g1, 4, 0,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_small(&g2, 4, 0,
                 1, 0, 1, 2, 2, 3, 3, 0, -1);

    /* a graph is always same as itself */
    err = igraph_is_same_graph(&g1, &g1, &res);
    IGRAPH_ASSERT(!err);
    IGRAPH_ASSERT(res);

    /* undirected graphs should be the same no matter
     * the direction of the edges (one is swapped in g2 */
    err = igraph_is_same_graph(&g1, &g2, &res);
    IGRAPH_ASSERT(!err);
    IGRAPH_ASSERT(res);    

    /* end of undirected */
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    /* directed graphs */
    igraph_small(&g1, 4, 1,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_small(&g2, 4, 1,
                 1, 0, 1, 2, 2, 3, 3, 0, -1);

    /* directed graphs should not be the same if an
     * edge has the opposite direction */
    err = igraph_is_same_graph(&g1, &g2, &res);
    IGRAPH_ASSERT(!err);
    IGRAPH_ASSERT(!res);

    igraph_destroy(&g2);

    /* change order of edges, they should be reordered by graph->ii */
    igraph_small(&g2, 4, 1,
                 1, 2, 0, 1, 2, 3, 3, 0, -1);
    err = igraph_is_same_graph(&g1, &g2, &res);
    IGRAPH_ASSERT(!err);
    IGRAPH_ASSERT(res);

    /* end of directed */
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    /* undirected vs directed */
    igraph_small(&g1, 4, 0,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    igraph_small(&g2, 4, 1,
                 0, 1, 1, 2, 2, 3, 3, 0, -1);
    err = igraph_is_same_graph(&g1, &g2, &res);
    IGRAPH_ASSERT(!err);
    IGRAPH_ASSERT(!res);

    /* end of undirected vs directed */
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    VERIFY_FINALLY_STACK();

    return 0;
}
