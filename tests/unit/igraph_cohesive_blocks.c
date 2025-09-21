/* igraph library.
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

void doit(igraph_t *g) {

    igraph_vector_int_list_t blocks;
    igraph_vector_int_t cohesion;
    igraph_vector_int_t parent;
    igraph_t block_tree;
    igraph_int_t i;

    igraph_vector_int_list_init(&blocks, 0);
    igraph_vector_int_init(&cohesion, 0);
    igraph_vector_int_init(&parent, 0);

    igraph_cohesive_blocks(g, &blocks, &cohesion, &parent,
                           &block_tree);

    printf("Blocks:\n");
    for (i = 0; i < igraph_vector_int_list_size(&blocks); i++) {
        igraph_vector_int_t *sg = igraph_vector_int_list_get_ptr(&blocks, i);
        printf("  ");
        igraph_vector_int_print(sg);
    }
    printf("Cohesion:\n  ");
    igraph_vector_int_print(&cohesion);
    printf("Parents:\n  ");
    igraph_vector_int_print(&parent);
    printf("Block graph:\n");
    igraph_write_graph_edgelist(&block_tree, stdout);

    igraph_vector_int_list_destroy(&blocks);
    igraph_vector_int_destroy(&cohesion);
    igraph_vector_int_destroy(&parent);
    igraph_destroy(&block_tree);

}

int main(void) {

    igraph_t g;

    printf("The graph from the Moody-White paper:\n");
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

    doit(&g);
    igraph_destroy(&g);
    printf("--\n");

    printf("Tricky graph, where the separators themselves form a block:\n");

    igraph_small(&g, 8, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 5, 1, 2, 1, 4, 1, 5, 1, 6, 2, 3, 2, 5, 2, 6, 2, 7,
                 3, 6, 3, 7, 4, 5, 5, 6, 6, 7,
                 -1);

    doit(&g);
    igraph_destroy(&g);
    printf("--\n");

    printf("The science camp graph from http://intersci.ss.uci.edu/wiki/index.php/Cohesive_blocking\n");
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

    doit(&g);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
