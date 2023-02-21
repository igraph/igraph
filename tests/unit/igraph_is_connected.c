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

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_bool_t conn;

    /* Null graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(! conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Singleton graph */
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(conn);

    igraph_destroy(&graph);

    /* Two isolated vertices, one with a self-loop */
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,0, -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(! conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Two isolated vertices, three self-loops */
    igraph_small(&graph, 2, IGRAPH_DIRECTED,
                 0,0, 0,0, 1,1, -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(! conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Weakly connected directed */
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 2,0, 1,2, 3,2,
                 -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(! conn);

    igraph_destroy(&graph);

    /* Directed cycle */
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 2,0, 1,3, 3,2,
                 -1);

    igraph_is_connected(&graph, &conn, IGRAPH_WEAK);
    IGRAPH_ASSERT(conn);

    igraph_is_connected(&graph, &conn, IGRAPH_STRONG);
    IGRAPH_ASSERT(conn);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
