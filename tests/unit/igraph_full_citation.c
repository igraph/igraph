/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <assert.h>
#include "test_utilities.inc"

int main() {
    igraph_t g, g_test;
    igraph_bool_t same;
    long int n_vertices = 4;

    /*    Undirected, should be a full graph    */
    assert(igraph_full_citation(&g, n_vertices, 0 /*undirected*/) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 4, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3, -1);
    assert(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    assert(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*    Directed, only edges from i->j if i > j    */
    assert(igraph_full_citation(&g, n_vertices, 1 /*directed*/) == IGRAPH_SUCCESS);
    igraph_small(&g_test, 4, IGRAPH_DIRECTED, 1, 0, 2, 0, 3, 0, 2, 1, 3, 1, 3, 2, -1);
    assert(igraph_is_same_graph(&g, &g_test, &same) == IGRAPH_SUCCESS);
    assert(same);
    igraph_destroy(&g);
    igraph_destroy(&g_test);

    /*    Directed, 1 vertex, should be edgeless    */
    assert(igraph_full_citation(&g, 1 /*n_vertices*/, 1 /*directed*/) == IGRAPH_SUCCESS);
    assert(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /*    Directed, 0 vertices, empty graph    */
    assert(igraph_full_citation(&g, 0 /*n_vertices*/, 1 /*directed*/) == IGRAPH_SUCCESS);
    assert(igraph_vcount(&g) == 0);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;

}
