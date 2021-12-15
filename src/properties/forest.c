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

/* count the number of vertices reachable from the root */
static int igraph_i_is_forest_visitor(igraph_integer_t root,igraph_vector_t *visited, const igraph_adjlist_t *al, igraph_integer_t *visited_count,igraph_integer_t *res,  igraph_neimode_t mode){

    igraph_stack_int_t stack;
    igraph_integer_t i;
    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);
    /* push the root into the stack */
    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));

    while (! igraph_stack_int_empty(&stack)) {
        igraph_integer_t u;
        igraph_vector_int_t *neighbors;
        igraph_integer_t ncount;
        u = igraph_stack_int_pop(&stack);

        /* Take a vertex from stack and check if it is already visited
         * if yes, then graph is not a forest
         * else add it to the visited vector
         */
        if (IGRAPH_LIKELY(! VECTOR(*visited)[u])) {
            VECTOR(*visited)[u] = 1;
            *visited_count += 1;
        }
        else{
            *res=0;
            break;
        }
        /* register all its neighbours (except its parent) for future processing */
        neighbors = igraph_adjlist_get(al, u);
        ncount = igraph_vector_int_size(neighbors);

        for (i = 0; i < ncount; ++i) {
            igraph_integer_t v = VECTOR(*neighbors)[i];
            if(mode==IGRAPH_ALL){
                if (IGRAPH_LIKELY(!VECTOR(*visited)[v])){
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, v));
                }
                /*to check for loop in undirected graph*/
                else if(v==u){
                    *res=0;
                    break;
                }
            }
            else{
                IGRAPH_CHECK(igraph_stack_int_push(&stack, v));
            }
        }
    }
    igraph_stack_int_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_is_forest
 * \brief Decides whether the graph is a forest.
 *
 * An undirected graph is a forest if it has no cycles.
 * </para><para>
 *
 * In the directed case, a possible additional requirement is that two
 * different trees should not have common vertex. Also, edges in each
 * tree are oriented away from the root (out-tree or arborescence) or all edges
 * are oriented towards the root (in-tree or anti-arborescence).
 * This test can be controlled using the \p mode parameter.
 * </para><para>
 *
 * By convention, the null graph (i.e. the graph with no vertices) is considered to be a forest.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored
 *        here.
 * \param roots A vector in which the root nodes will be stored. When \p mode
 *        is \c IGRAPH_ALL or the graph is undirected, any 1 vertex from each
 *        component can be the root. When \p mode is \c IGRAPH_OUT
 *        or \c IGRAPH_IN, all the vertices with zero in- or out-degree,
 *        respectively are considered as root nodes.
 * \param mode For a directed graph this specifies whether to test for an
 *        out-tree, an in-tree or ignore edge directions. The respective
 *        possible values are:
 *        \c IGRAPH_OUT, \c IGRAPH_IN, \c IGRAPH_ALL. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: At most O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 */
int igraph_is_forest(const igraph_t *graph,igraph_bool_t *res, igraph_vector_t *roots , igraph_neimode_t mode) {
    igraph_adjlist_t al;
    igraph_vector_t visited;
    igraph_integer_t visited_count=0;
    igraph_integer_t vcount, ecount;

    vcount = igraph_vcount(graph);
    ecount = igraph_ecount(graph);

    /*By convention, a zero-vertex graph will be considered a forest.*/
    if (vcount == 0) {
        *res = 1;
        return IGRAPH_SUCCESS;
    }
    // A forest can have maximum vcount-1 edges.
    if (ecount > vcount - 1) {
        *res = 0;
        return IGRAPH_SUCCESS;
    }
    /*A single-vertex graph is a forest, provided it has no edges (checked in the previous if (..)) */
    if (vcount == 1) {
        *res = 1;
        igraph_vector_push_back(roots,0);
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

    IGRAPH_CHECK(igraph_vector_resize(&visited,vcount));
    igraph_integer_t i;
    for(i =0; i<vcount;++i){
        VECTOR(visited)[i]=0;
    }

    /* The main algorithm:
     * Undirected Graph:- We add each unvisited vertex to the roots vector, and
     * mark all other vertices that are reachable from it as visited.
     *
     * Directed Graph:- For each tree, the root is the node with no
     * incoming/outgoing connections, depending on 'mode'. We add each vertex
     * with zero degree to the roots vector and mark all other vertices that are
     * reachable from it as visited.
     *
     * If all the vertices are visited exactly once, then the graph is a forest.
     */

    switch (mode) {
    case IGRAPH_ALL:{

        for(i =0; i<vcount;++i){
            if (!VECTOR(visited)[i]) {
                igraph_vector_push_back(roots,i);
                IGRAPH_CHECK(igraph_i_is_forest_visitor(i,&visited, &al, &visited_count,res,mode));
            }
        }
        break;
    }
    case IGRAPH_IN:
    case IGRAPH_OUT: {
        igraph_vector_t degree;

        IGRAPH_CHECK(igraph_vector_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_destroy, &degree);
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), mode == IGRAPH_IN ? IGRAPH_OUT : IGRAPH_IN, /* loops = */ 1));

        for (i = 0; i < vcount; ++i){
            if (VECTOR(degree)[i] == 0) {
                igraph_vector_push_back(roots,i);
                IGRAPH_CHECK(igraph_i_is_forest_visitor(i,&visited, &al, &visited_count,res,mode));
            }
        }

        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
        break;
    }
    default:
        IGRAPH_ERROR("Invalid mode", IGRAPH_EINVMODE);
    }
    if(*res){
        *res= visited_count==vcount;
    }
    /*If the graph is not a forest then the root vector will be empty*/
    if(!*res){
        IGRAPH_CHECK(igraph_vector_resize(roots,0));
    }

    igraph_vector_destroy(&visited);
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
