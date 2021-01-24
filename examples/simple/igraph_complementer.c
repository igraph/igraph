/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

    igraph_t g1, g2;

    /* complementer of the empty graph */
    igraph_empty(&g1, 5, IGRAPH_DIRECTED);
    igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
    igraph_write_graph_edgelist(&g2, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /* the same without loops */
    igraph_empty(&g1, 5, IGRAPH_DIRECTED);
    igraph_complementer(&g2, &g1, IGRAPH_NO_LOOPS);
    igraph_write_graph_edgelist(&g2, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /* complementer of the full graph */
    igraph_full(&g1, 5, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
    if (igraph_ecount(&g2) != 0) {
        return 1;
    }
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /* complementer of the full graph, results loops only */
    igraph_full(&g1, 5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
    igraph_write_graph_edgelist(&g2, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /**************
     * undirected *
     *************/

    /* complementer of the empty graph */
    igraph_empty(&g1, 5, IGRAPH_UNDIRECTED);
    igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
    igraph_write_graph_edgelist(&g2, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /* the same without loops */
    igraph_empty(&g1, 5, IGRAPH_UNDIRECTED);
    igraph_complementer(&g2, &g1, IGRAPH_NO_LOOPS);
    igraph_write_graph_edgelist(&g2, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /* complementer of the full graph */
    igraph_full(&g1, 5, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
    if (igraph_ecount(&g2) != 0) {
        return 1;
    }
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    printf("---\n");

    /* complementer of the full graph, results loops only */
    igraph_full(&g1, 5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_complementer(&g2, &g1, IGRAPH_LOOPS);
    igraph_write_graph_edgelist(&g2, stdout);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    return 0;
}
