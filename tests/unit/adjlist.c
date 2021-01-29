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

int test_simple_trees() {
    igraph_t g, g2;
    igraph_adjlist_t adjlist;
    igraph_bool_t iso;

    /* Directed, out */
    igraph_tree(&g, 42, 3, IGRAPH_TREE_OUT);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_OUT, /*duplicate=*/ 0);
    igraph_isomorphic(&g, &g2, &iso);
    if (!iso) {
        return 1;
    }
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    /* Directed, in */
    igraph_tree(&g, 42, 3, IGRAPH_TREE_OUT);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_IN, /*duplicate=*/ 0);
    igraph_isomorphic(&g, &g2, &iso);
    if (!iso) {
        return 1;
    }
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    /* Undirected */
    igraph_tree(&g, 42, 3, IGRAPH_TREE_UNDIRECTED);
    igraph_adjlist_init(&g, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE);
    igraph_adjlist(&g2, &adjlist, IGRAPH_ALL, /*duplicate=*/ 1);
    igraph_isomorphic(&g, &g2, &iso);
    if (!iso) {
        return 1;
    }
    igraph_adjlist_destroy(&adjlist);
    igraph_destroy(&g2);
    igraph_destroy(&g);

    return 0;
}

int main() {
    int retval;

    RUN_TEST(test_simple_trees);

    return 0;
}
