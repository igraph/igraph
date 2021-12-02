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

int main() {
    
    igraph_t graph;
    igraph_bool_t res;
    igraph_vector_int_t v;
    igraph_vector_int_init_int(&v, 3, 3, 4, 5);

    /* Create a directed symmetric tree with 2 levels - 
       3 children in first and 4 children in second level,
       5 children in third level
       with edges pointing towards the root. */
    igraph_symmetric_tree(&graph, &v, IGRAPH_TREE_IN);

    igraph_is_tree(&graph, &res, NULL, IGRAPH_IN);
    printf("Is it an in-tree? %s\n", res ? "Yes" : "No");

    igraph_is_tree(&graph, &res, NULL, IGRAPH_OUT);
    printf("Is it an out-tree? %s\n", res ? "Yes" : "No");

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&v);

    return 0;
}