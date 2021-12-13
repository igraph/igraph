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

#include "igraph_structural.h"

#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_stack.h"

static int igraph_i_is_forest_visitor(igraph_integer_t root,const igraph_adjlist_t *al, igraph_vector_t *visited,  igraph_integer_t *visited_count,igraph_integer_t *res,  igraph_neimode_t mode){
    igraph_stack_int_t stack;
    long i;
    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);
    /* push the root into the stack */
    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));

    // implement function

    igraph_stack_int_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(1);
    
    return IGRAPH_SUCCESS;
}


int igraph_is_forest(const igraph_t *graph,igraph_bool_t *res, igraph_vector_t *roots , igraph_neimode_t mode) {
    igraph_adjlist_t al;
    igraph_vector_t visited;
    igraph_integer_t visited_count=0;
    igraph_integer_t vcount, ecount;
    
    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    /* A forest can have maximum vcount-1 edges. */
    if (ecount > vcount - 1) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }
    /* The single-vertex graph is a forest, provided it has no edges (checked in the previous if (..)) */
    if (vcount == 1) {
        *res = 1;
        return IGRAPH_SUCCESS;
    }
    /* Ignore mode for undirected graphs. */
    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }
    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    
    *res = 1; /* assume success */
    IGRAPH_CHECK(igraph_vector_init(&visited, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &visited);
    
    switch (mode) {
    case IGRAPH_ALL:
        igraph_integer_t i;
        
        //for undirected graph

        break;

    case IGRAPH_IN:
    case IGRAPH_OUT: {
        igraph_vector_t degree;

        IGRAPH_CHECK(igraph_vector_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_destroy, &degree);

        // for directed graph 

        break;
    }
    default:
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVMODE);
    }
    
    igraph_vector_destroy(&visited);
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
