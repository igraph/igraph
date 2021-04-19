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

int main() {

    igraph_t g;

    /* Multiple edges */

    igraph_small(&g, 0, IGRAPH_DIRECTED, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, -1);
    igraph_simplify(&g, 1, 1, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, -1);
    igraph_simplify(&g, 1, 1, /*edge_comb=*/ 0);
    if (igraph_ecount(&g) != 1) {
        return 1;
    }
    igraph_destroy(&g);

    /* Loop edges*/

    igraph_small(&g, 0, IGRAPH_DIRECTED, 0, 0, 1, 1, 2, 2, 1, 2, -1);
    igraph_simplify(&g, 1, 1, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 0, 0, 1, 1, 2, 2, 1, 2, -1);
    igraph_simplify(&g, 1, 1, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    /* Loop & multiple edges */

    igraph_small(&g, 0, IGRAPH_DIRECTED, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, -1);
    igraph_simplify(&g, 1 /* multiple */, 0 /* loop */, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, -1);
    igraph_simplify(&g, 1 /* multiple */, 0 /* loop */, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_DIRECTED, 2, 2, 2, 2, 2, 2, 3, 2, -1);
    igraph_simplify(&g, 0 /* multiple */, 1 /* loop */, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 3, 3, 3, 3, 3, 4, -1);
    igraph_simplify(&g, 0 /* multiple */, 1 /* loop */, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_DIRECTED, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, -1);
    igraph_simplify(&g, 1, 1, /*edge_comb=*/ 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 2, 3, 2, -1);
    igraph_simplify(&g, 1, 1, /*edge_comb=*/ 0);
    if (igraph_ecount(&g) != 1) {
        return 2;
    }
    igraph_destroy(&g);

    return 0;
}
