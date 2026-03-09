/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

void check_simple(const igraph_t *graph,
                  igraph_bool_t expected_directed,
                  igraph_bool_t expected_undirected) {

    igraph_t ug;
    igraph_bool_t simple;
    igraph_adjlist_t al;

    igraph_is_simple(graph, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_directed) );

    igraph_is_simple(graph, &simple, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_undirected) );

    igraph_copy(&ug, graph);
    igraph_to_undirected(&ug, IGRAPH_TO_UNDIRECTED_EACH, NULL);
    igraph_is_simple(&ug, &simple, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_undirected) );
    igraph_destroy(&ug);

    /* Reset the cache and try again. This clears cache entries that may have
     * been set during graph construction. */

    igraph_invalidate_cache(graph);

    igraph_is_simple(graph, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_directed) );

    igraph_invalidate_cache(graph);

    igraph_is_simple(graph, &simple, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_undirected) );

    /* Reset the cache, run adjlist constructor (which does set the cache)
     * and try again. This is meant to verify that adjlist construcors,
     * which are used particularly frequently, do not corrupt the cache. */

    /* Undirected adjlist for directed graph. */
    igraph_invalidate_cache(graph);
    igraph_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    igraph_adjlist_destroy(&al);

    igraph_is_simple(graph, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_directed) );

    /* Directed adjlist for directed graph. */
    igraph_invalidate_cache(graph);
    igraph_adjlist_init(graph, &al, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    igraph_adjlist_destroy(&al);

    igraph_is_simple(graph, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_directed) );

    /* Undirected adjlist for directed graph. */
    igraph_invalidate_cache(graph);
    igraph_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    igraph_adjlist_destroy(&al);

    igraph_is_simple(graph, &simple, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_undirected) );

    /* Directed adjlist for directed graph. */
    igraph_invalidate_cache(graph);
    igraph_adjlist_init(graph, &al, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE);
    igraph_adjlist_destroy(&al);

    igraph_is_simple(graph, &simple, IGRAPH_UNDIRECTED);
    IGRAPH_ASSERT( (!simple) == (!expected_undirected) );
}

int main(void) {
    igraph_t graph;

    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    check_simple(&graph, true, true);
    igraph_destroy(&graph);

    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    check_simple(&graph, true, true);
    igraph_destroy(&graph);

    igraph_empty(&graph, 2, IGRAPH_DIRECTED);
    check_simple(&graph, true, true);
    igraph_destroy(&graph);

    igraph_full(&graph, 1, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    check_simple(&graph, false, false);
    igraph_destroy(&graph);

    igraph_full(&graph, 2, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    check_simple(&graph, true, false);
    igraph_destroy(&graph);

    igraph_full(&graph, 2, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    check_simple(&graph, false, false);
    igraph_destroy(&graph);

    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 3,0, 0,3,
                 -1);
    check_simple(&graph, true, false);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
