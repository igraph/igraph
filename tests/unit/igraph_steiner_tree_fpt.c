 
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

int main() {
    igraph_t g_empty, g_empty_dir, g_lm;
    
    igraph_vector_t steiner_terminals;
    igraph_vector_t weights_empty, weights_lm;
    igraph_integer_t dreyfus_wagner_out;
   

    igraph_vector_init(&weights_empty, 0);
    
    igraph_vector_init(&steiner_terminals,0,2,4);


    igraph_vector_init_int(&weights_lm ,3, 1, 7, 5, 1, 2, 7);
    
    igraph_small(&g_empty, 0, 0, -1);
 
    igraph_small(&g_lm, 5,IGRAPH_UNDIRECTED,0,1, 0,2, 1,2, 1,3, 1,4, 2,3, 3,4 -1);

  
    printf("No vertices, not directed:\n");
    IGRAPH_ASSERT(igraph_steiner_dreyfus_wagner(&g_empty, &steiner_terminals,IGRAPH_ALL,&weights_empty ) == IGRAPH_SUCCESS);
    
    
    
    printf("Un-Directed graph with loops and multi-edges, select none:\n");
    dreyfus_wagner_out = igraph_steiner_dreyfus_wagner(&g_lm, &steiner_terminals,IGRAPH_ALL,&weights_lm);
  
    printf(dreyfus_wagner_out);
    
    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);
    igraph_vector_destroy(&weights_empty);
    igraph_vector_destroy(&weights_lm);
  

    VERIFY_FINALLY_STACK();
    return 0;
}

   
