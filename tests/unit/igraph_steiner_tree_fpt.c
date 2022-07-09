
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
#include "test_utilities.h"

int main()
{
   // printf("Starting Tests\n");
    igraph_t g_empty, g_lm;

    igraph_vector_int_t steiner_terminals,steiner_terminals_null;
    igraph_vector_t weights_empty, weights_lm;
    
    igraph_vector_init(&weights_empty, 0);

    igraph_vector_int_init(&steiner_terminals, 4);
    igraph_vector_int_init(&steiner_terminals_null, 0);

    VECTOR(steiner_terminals)[0] = 0;
    VECTOR(steiner_terminals)[1] = 1;
    VECTOR(steiner_terminals)[2] = 2;
    VECTOR(steiner_terminals)[3] = 3;

    igraph_vector_init(&weights_lm, 21);
    
    VECTOR(weights_lm)[0] = 2;
    VECTOR(weights_lm)[1] = 2;
    VECTOR(weights_lm)[2] = 2;
    VECTOR(weights_lm)[3] = 1;
    VECTOR(weights_lm)[4] = 1;
    VECTOR(weights_lm)[5] = 2;
    
    VECTOR(weights_lm)[6] = 2;
    VECTOR(weights_lm)[7] = 2;
    VECTOR(weights_lm)[8] = 2;
    VECTOR(weights_lm)[9] = 1;
    VECTOR(weights_lm)[10] = 2;

    VECTOR(weights_lm)[11] = 2;
    VECTOR(weights_lm)[12] = 2;
    VECTOR(weights_lm)[13] = 2;
    VECTOR(weights_lm)[14] = 1;
    
    VECTOR(weights_lm)[15] = 1;
    VECTOR(weights_lm)[16] = 2;
    VECTOR(weights_lm)[17] = 1;

    VECTOR(weights_lm)[18] = 2;
    VECTOR(weights_lm)[19] = 1;

    VECTOR(weights_lm)[20] = 1;

    igraph_empty(&g_empty, 0, 0);

    igraph_small(&g_lm, 7, IGRAPH_UNDIRECTED, 
                                                0, 1,
                                                0, 2,
                                                0, 3, 
                                                0, 4, 
                                                0, 5, 
                                                0, 6, 
                                                
                                                1, 2, 
                                                1, 3, 
                                                1, 4, 
                                                1, 5, 
                                                1, 6,

                                                2, 3, 
                                                2, 4, 
                                                2, 5, 
                                                2, 6, 
                                                
                                                3, 4, 
                                                3, 5, 
                                                3, 6, 
                                                
                                                4, 5, 
                                                4, 6, 
                                                
                                                5, 6, 
                                                -1);
    
    printf("No vertices, not directed:\n");
    igraph_real_t val1,val2;
    IGRAPH_ASSERT(igraph_steiner_dreyfus_wagner(&g_empty,&steiner_terminals_null, IGRAPH_ALL, &weights_empty,&val1) == IGRAPH_FAILURE);
    printf("%.2f\n",val1);
    IGRAPH_ASSERT(val1 == 0);
    printf("Un-Directed graph with loops and multi-edges, select none:\n");
    IGRAPH_ASSERT(igraph_steiner_dreyfus_wagner(&g_lm,&steiner_terminals, IGRAPH_ALL, &weights_lm,&val2) == IGRAPH_SUCCESS);
    printf("%.2f\n",val2);
    IGRAPH_ASSERT(val2 == 5);
    igraph_destroy(&g_empty);
    igraph_destroy(&g_lm);
    igraph_vector_destroy(&weights_empty);
    igraph_vector_destroy(&weights_lm);

    igraph_vector_int_destroy(&steiner_terminals);
    igraph_vector_int_destroy(&steiner_terminals_null);

   

    VERIFY_FINALLY_STACK();
    
    return 0;
}
