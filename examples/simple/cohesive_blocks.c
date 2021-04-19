/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int doit(igraph_t *g) {

    igraph_vector_ptr_t blocks;
    igraph_vector_t cohesion;
    igraph_vector_t parent;
    igraph_t block_tree;
    long int i;

    igraph_vector_ptr_init(&blocks, 0);
    igraph_vector_init(&cohesion, 0);
    igraph_vector_init(&parent, 0);

    igraph_cohesive_blocks(g, &blocks, &cohesion, &parent,
                           &block_tree);

    printf("Blocks:\n");
    for (i = 0; i < igraph_vector_ptr_size(&blocks); i++) {
        igraph_vector_t *sg = VECTOR(blocks)[i];
        printf("  ");
        igraph_vector_print(sg);
        igraph_vector_destroy(sg);
        igraph_free(sg);
    }
    printf("Cohesion:\n  ");
    igraph_vector_print(&cohesion);
    printf("Parents:\n  ");
    igraph_vector_print(&parent);
    printf("Block graph:\n");
    igraph_write_graph_edgelist(&block_tree, stdout);

    igraph_vector_ptr_destroy(&blocks);
    igraph_vector_destroy(&cohesion);
    igraph_vector_destroy(&parent);
    igraph_destroy(&block_tree);

    return 0;
}

int main() {

    igraph_t g;
    int ret;

    /* --------------------------------------------------------*/
    /* The graph from the Moody-White paper                    */

    igraph_small(&g, 23, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5,
                 1, 2, 1, 3, 1, 4, 1, 6,
                 2, 3, 2, 5, 2, 6,
                 3, 4, 3, 5, 3, 6,
                 4, 5, 4, 6, 4, 20,
                 5, 6,
                 6, 7, 6, 10, 6, 13, 6, 18,
                 7, 8, 7, 10, 7, 13,
                 8, 9,
                 9, 11, 9, 12,
                 10, 11, 10, 13,
                 11, 15,
                 12, 15,
                 13, 14,
                 14, 15,
                 16, 17, 16, 18, 16, 19,
                 17, 19, 17, 20,
                 18, 19, 18, 21, 18, 22,
                 19, 20,
                 20, 21, 20, 22,
                 21, 22,
                 -1);

    if ( (ret = doit(&g)) ) {
        return ret;
    }
    igraph_destroy(&g);
    printf("--\n");

    /* --------------------------------------------------------*/
    /* A tricky graph, where the separators themselves         */
    /* form a block. But recently we don't include this        */
    /* block in the results.                                   */

    igraph_small(&g, 8, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 5, 1, 2, 1, 4, 1, 5, 1, 6, 2, 3, 2, 5, 2, 6, 2, 7,
                 3, 6, 3, 7, 4, 5, 5, 6, 6, 7,
                 -1);

    if ( (ret = doit(&g)) ) {
        return ret;
    }
    igraph_destroy(&g);
    printf("--\n");

    /* --------------------------------------------------------*/
    /* The science camp graph from http://intersci.ss.uci.edu/ */
    /* wiki/index.php/Cohesive_blocking                        */

    igraph_small(&g, 18, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3,
                 1, 2, 1, 3, 1, 16, 1, 17,
                 2, 3,
                 3, 17,
                 4, 5, 4, 6, 4, 7, 4, 8,
                 5, 6, 5, 7,
                 6, 7, 6, 8,
                 7, 8, 7, 16,
                 8, 9, 8, 10,
                 9, 11, 9, 12, 9, 13, 9, 14,
                 10, 11, 10, 12, 10, 13,
                 11, 14,
                 12, 13, 12, 14, 12, 15,
                 15, 16, 15, 17,
                 16, 17,
                 -1);

    if ( (ret = doit(&g)) ) {
        return ret;
    }
    igraph_destroy(&g);
    printf("--\n");

    /* --------------------------------------------------------*/
    /* Zachary karate-club                                     */

    igraph_small(&g, 34, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 10, 0, 11, 0, 12, 0, 13,
                 0, 17, 0, 19, 0, 21, 0, 31,
                 1, 2, 1, 3, 1, 7, 1, 13, 1, 17, 1, 19, 1, 21, 1, 30,
                 2, 3, 2, 7, 2, 27, 2, 28, 2, 32, 2, 9, 2, 8, 2, 13,
                 3, 7, 3, 12, 3, 13,
                 4, 6, 4, 10,
                 5, 6, 5, 10, 5, 16,
                 6, 16,
                 8, 30, 8, 32, 8, 33,
                 9, 33,
                 13, 33,
                 14, 32, 14, 33,
                 15, 32, 15, 33,
                 18, 32, 18, 33,
                 19, 33,
                 20, 32, 20, 33,
                 22, 32, 22, 33,
                 23, 25, 23, 27, 23, 32, 23, 33, 23, 29,
                 24, 25, 24, 27, 24, 31,
                 25, 31,
                 26, 29, 26, 33,
                 27, 33,
                 28, 31, 28, 33,
                 29, 32, 29, 33,
                 30, 32, 30, 33,
                 31, 32, 31, 33,
                 32, 33,
                 -1);

    if ( (ret = doit(&g)) ) {
        return ret;
    }
    igraph_destroy(&g);
    printf("--\n");

    return 0;
}
