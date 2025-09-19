/*
   igraph library.
   Copyright (C) 2009-2021  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t vids, layers, parents;

    igraph_vector_int_init(&vids, 0);
    igraph_vector_int_init(&layers, 0);
    igraph_vector_int_init(&parents, 0);

    /* Test a ring graph */
    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 0);
    igraph_bfs_simple(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);
    print_vector_int(&vids);
    print_vector_int(&layers);
    print_vector_int(&parents);
    igraph_destroy(&g);

    /* Test a tree graph */
    igraph_kary_tree(&g, 20, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_bfs_simple(&g, 0, IGRAPH_ALL, &vids, &layers, &parents);
    print_vector_int(&vids);
    print_vector_int(&layers);
    print_vector_int(&parents);
    igraph_destroy(&g);

    /* Test th same graph with all arguments as nulls to see if we tolerate that */
    igraph_kary_tree(&g, 20, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_bfs_simple(&g, 0, IGRAPH_ALL, 0, 0, 0);
    igraph_destroy(&g);

    /* Also test the case when 'layers' is not null and 'vids' is null to ensure
     * that we don't need 'vids' in the internal implementation to populate
     * 'layers' */
    igraph_kary_tree(&g, 20, 2, IGRAPH_TREE_UNDIRECTED);
    igraph_bfs_simple(&g, 0, IGRAPH_ALL, 0, &layers, 0);
    print_vector_int(&layers);
    igraph_destroy(&g);

    /* Test directed graph where not all nodes are reachable */
    igraph_kary_tree(&g, 20, 2, IGRAPH_TREE_OUT);
    igraph_bfs_simple(&g, 7, IGRAPH_OUT, 0, &layers, &parents);
    print_vector_int(&layers);
    print_vector_int(&parents);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&vids);
    igraph_vector_int_destroy(&layers);
    igraph_vector_int_destroy(&parents);

    VERIFY_FINALLY_STACK();

    return 0;
}
