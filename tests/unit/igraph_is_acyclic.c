/*
   igraph library.
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
    igraph_bool_t acyclic;

    /* Null graph directed */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Null graph undirected */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Singleton graph directed */
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Singleton graph undirected */
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Directed cyclic */
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
        0,1, 2,0, 1,3, 3,2,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(!acyclic);
    igraph_destroy(&graph);

    /* Directed acyclic */
    igraph_small(&graph, 3, IGRAPH_DIRECTED,
        0,1, 2,0, 1,3,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Undirected cyclic */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
        0,1, 2,0, 1,3, 3,2,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(! acyclic);
    igraph_destroy(&graph);

    /* Undirected acyclic */
    igraph_small(&graph, 3, IGRAPH_UNDIRECTED,
        0,1, 2,0, 1,3,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    /* Self loop directed */
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
        0,1, 1,3, 3,2, 2,2,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(! acyclic);
    igraph_destroy(&graph);

    /* Self loop undirected */
    igraph_small(&graph, 4, IGRAPH_UNDIRECTED,
        0,1, 2,0, 1,3, 3,3,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(! acyclic);
    igraph_destroy(&graph);

    /* Directed acyclic graph which would be cyclic if undirected */
    igraph_small(&graph, 3, IGRAPH_DIRECTED,
        0,1, 1,2, 0,2,
        -1);
    igraph_is_acyclic(&graph, &acyclic);
    IGRAPH_ASSERT(acyclic);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;

}
