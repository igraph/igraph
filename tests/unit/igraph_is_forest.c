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

void check_output(const igraph_t *graph,igraph_bool_t *res,igraph_neimode_t mode){
    igraph_bool_t result;
    igraph_vector_t roots;
    igraph_vector_init(&roots, 0);
    IGRAPH_FINALLY(igraph_vector_destroy, &roots);
    igraph_is_forest(graph,&result,&roots, mode);
    IGRAPH_ASSERT(*res==result);
    if(result)
    igraph_vector_print(&roots);
    igraph_vector_destroy(&roots);
}

int main() {
    igraph_t graph;
    igraph_bool_t res;
    igraph_neimode_t mode;

    /*Empty Graph*/
    mode=IGRAPH_ALL;
    res=1;
    igraph_small(&graph,2,0,-1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    /*Graph with 0 edges*/
    mode=IGRAPH_ALL;
    res=1;
    igraph_small(&graph,5,0,-1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    /*Undirected Graph*/
    mode=IGRAPH_ALL;
    res=1;
    igraph_small(&graph,5,0, 0,1, 1,2, 3,4, 3,5, -1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    /*Directed Graph out trees*/
    mode=IGRAPH_OUT;
    res=1;
    igraph_small(&graph,5,1, 0,1, 1,2, 3,4, 3,5, -1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    /*Directed Graph in trees*/
    mode=IGRAPH_OUT;
    res=1;
    igraph_small(&graph,5,1, 0,1, 1,2, 3,4, 3,5, -1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    /*Undirected Graph with cycle*/
    mode=IGRAPH_ALL;
    res=0;
    igraph_small(&graph,5,0, 0,1, 1,2, 3,4, 3,5, 4,5 -1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    mode=IGRAPH_ALL;
    res=0;
    igraph_small(&graph,4,0, 0,1, 1,2, 2,3, 3,3, -1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    /*Directed Graph*/
    mode=IGRAPH_OUT;
    res=0;
    igraph_small(&graph,5,1, 0,1, 1,2, 3,4, 3,5, 4,2, -1);
    check_output(&graph,&res,mode);
    igraph_destroy(&graph);

    //VERIFY_FINALLY_STACK();

  return 0;
}
