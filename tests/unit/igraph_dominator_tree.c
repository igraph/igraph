/*
   IGraph library.
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

#include <igraph.h>
#include "test_utilities.h"

int main(void) {
    igraph_t g, domtree;
    igraph_vector_int_t dom;
    igraph_vector_int_t leftout;

    igraph_vector_int_init(&dom, 0);
    igraph_small(&g, 13, IGRAPH_DIRECTED,
                 0, 1, 0, 7, 0, 10,
                 1, 2, 1, 5,
                 2, 3,
                 3, 4,
                 4, 3, 4, 0,
                 5, 3, 5, 6,
                 6, 3,
                 7, 8, 7, 10, 7, 11,
                 8, 9,
                 9, 4, 9, 8,
                 10, 11,
                 11, 12,
                 12, 9,
                 -1);

    /* Check NULL vector arguments */
    igraph_dominator_tree(&g, /*root=*/ 0, /*dom=*/ 0, /*domtree=*/ 0,
                          /*leftout=*/ 0, /*mode=*/ IGRAPH_OUT);

    /* Proper calculation */
    igraph_dominator_tree(&g, /*root=*/ 0, &dom, /*domtree=*/ 0,
                          /*leftout=*/ 0, /*mode=*/ IGRAPH_OUT);
    igraph_vector_int_print(&dom);

    /* Tree calculation */
    igraph_dominator_tree(&g, /*root=*/ 0, /*dom=*/ 0, /*domtree=*/ &domtree,
                          /*leftout=*/ 0, /*mode=*/ IGRAPH_OUT);
    igraph_write_graph_edgelist(&domtree, stdout);

    igraph_vector_int_destroy(&dom);
    igraph_destroy(&domtree);
    igraph_destroy(&g);

    /* -------------------------------------------------------------------*/

    igraph_vector_int_init(&dom, 0);
    igraph_small(&g, 13, IGRAPH_DIRECTED,
                 1, 0, 2, 0, 3, 0,
                 4, 1,
                 1, 2, 4, 2, 5, 2,
                 6, 3, 7, 3,
                 12, 4,
                 8, 5,
                 9, 6,
                 9, 7, 10, 7,
                 5, 8, 11, 8,
                 11, 9,
                 9, 10,
                 9, 11, 0, 11,
                 8, 12,
                 -1);

    /* Proper calculation */
    igraph_dominator_tree(&g, /*root=*/ 0, &dom, /*domtree=*/ 0,
                          /*leftout=*/ 0, /*mode=*/ IGRAPH_IN);
    igraph_vector_int_print(&dom);

    /* Tree calculation */
    igraph_dominator_tree(&g, /*root=*/ 0, /*dom=*/ 0, /*domtree=*/ &domtree,
                          /*leftout=*/ 0, /*mode=*/ IGRAPH_IN);
    igraph_write_graph_edgelist(&domtree, stdout);

    igraph_vector_int_destroy(&dom);
    igraph_destroy(&domtree);
    igraph_destroy(&g);

    /* -------------------------------------------------------------------*/

    igraph_vector_int_init(&dom, 0);
    igraph_vector_int_init(&leftout, 0);

    /* Check a graph with more components */
    igraph_small(&g, 20, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,
                 1, 4,
                 2, 1, 2, 4, 2, 8,
                 3, 9, 3, 10,
                 4, 15,
                 8, 11,
                 9, 12,
                 10, 12, 10, 13,
                 11, 8, 11, 14,
                 12, 14,
                 13, 12,
                 14, 12, 14, 0,
                 15, 11,
                 -1);

    igraph_dominator_tree(&g, /*root=*/ 0, &dom, &domtree,
                          &leftout, /*mode=*/ IGRAPH_OUT);
    igraph_vector_int_print(&dom);
    igraph_vector_int_print(&leftout);
    igraph_write_graph_edgelist(&domtree, stdout);

    igraph_vector_int_destroy(&dom);
    igraph_vector_int_destroy(&leftout);
    igraph_destroy(&domtree);
    igraph_destroy(&g);


    VERIFY_FINALLY_STACK();

    return 0;
}
