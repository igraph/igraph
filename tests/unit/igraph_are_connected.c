/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2021  The igraph development team <igraph@igraph.org>

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
    igraph_t g;
    igraph_bool_t connected;

    /* Complete graph. Any two distinct vertices are connected. */

    igraph_full(&g, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_are_connected(&g, 2, 7, &connected);
    IGRAPH_ASSERT(connected);

    igraph_are_connected(&g, 0, 0, &connected);
    IGRAPH_ASSERT(! connected);

    igraph_destroy(&g);

    /* Complete graph with self-loops. */

    igraph_full(&g, 10, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    igraph_are_connected(&g, 1, 7, &connected);
    IGRAPH_ASSERT(connected);

    igraph_are_connected(&g, 2, 2, &connected);
    IGRAPH_ASSERT(connected);

    igraph_destroy(&g);

    /* Graph with no edges. Any two distinct vertices are disconnected. */

    igraph_empty(&g, 10, IGRAPH_DIRECTED);
    igraph_are_connected(&g, 3, 6, &connected);
    IGRAPH_ASSERT(! connected);
    igraph_destroy(&g);

    /* Invalid vertex ID, expecting an error */

    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,0, 2,3, -1);
    igraph_set_error_handler(igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_are_connected(&g, 0, igraph_vcount(&g) + 2 /* vertex id out of range */, &connected) == IGRAPH_EINVVID);
    igraph_set_error_handler(igraph_error_handler_abort);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
