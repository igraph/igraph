/*
   igraph library.
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

int main(void) {
    igraph_t g;
    igraph_vector_int_list_t blocks;
    igraph_vector_int_t cohesion;
    igraph_vector_int_t parent;
    igraph_t block_tree;
    igraph_int_t i;

    /* Initialize the library. */
    igraph_setup();

    igraph_famous(&g, "zachary");
    igraph_vector_int_list_init(&blocks, 0);
    igraph_vector_int_init(&cohesion, 0);
    igraph_vector_int_init(&parent, 0);

    igraph_cohesive_blocks(&g, &blocks, &cohesion, &parent,
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

    igraph_destroy(&g);

    return 0;
}
