/*
   igraph library.
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

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_bool_t has_mutual;

    /* undirected null graph */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(! has_mutual);
    igraph_destroy(&graph);

    /* undirected edgeless graph */
    igraph_empty(&graph, 3, IGRAPH_UNDIRECTED);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(! has_mutual);
    igraph_destroy(&graph);

    /* directed null graph */
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(! has_mutual);
    igraph_destroy(&graph);

    /* directed edgeless graph */
    igraph_empty(&graph, 3, IGRAPH_DIRECTED);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(! has_mutual);
    igraph_destroy(&graph);

    /* undirected with edges */
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED, 0,1, -1);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(has_mutual);
    igraph_destroy(&graph);

    /* directed with no mutual */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,2, -1);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(! has_mutual);
    igraph_destroy(&graph);

    /* directed with mutual */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 1,0, -1);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(has_mutual);
    igraph_destroy(&graph);

    /* directed with loops, loops considered mutual */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 0,1, 1,1, 1,2, -1);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_LOOPS);
    IGRAPH_ASSERT(has_mutual);
    igraph_destroy(&graph);

    /* directed with loops, loops not considered mutual */
    igraph_small(&graph, 0, IGRAPH_DIRECTED, 0,1, 0,1, 1,1, 1,2, -1);
    igraph_has_mutual(&graph, &has_mutual, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(!has_mutual);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
