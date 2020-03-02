#include <stdio.h>
/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_structural.h"
#include "igraph_transitivity.h"
#include "igraph_paths.h"
#include "igraph_math.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_centrality.h"
#include "igraph_components.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_types_internal.h"
#include "igraph_dqueue.h"
#include "igraph_attributes.h"
#include "igraph_neighborhood.h"
#include "igraph_topology.h"
#include "igraph_qsort.h"
#include "config.h"
#include "structural_properties_internal.h"

#include <assert.h>

// solution adapted from https://www.geeksforgeeks.org/eulerian-path-and-circuit/
// The function returns one of the following values
// 0 -> if graph is not Eulerian
// 1 -> if graph has an euler circuit
// 2 -> if graph has an euler path

// only considering undirected graphs at the moment
int is_eulerian(igraph_t *graph) {

    igraph_bool_t res;
    igraph_is_connected(graph, &res, IGRAPH_WEAK);
    if (res == 0) return 0;

    igraph_inclist_t il;
    igraph_inclist_init(graph, &il, IGRAPH_ALL);
    igraph_vector_int_t *incedges;

    int odd = 0;
    for (int i = 0; i < igraph_vcount(graph); i++) {
        incedges = igraph_inclist_get(&il, i);
        if ( igraph_vector_int_size(incedges) % 2 == 1) odd++;
    }

    if (odd > 2) return 0;

    // if pdd count is 2, then semi-eulerian
    // if odd count is 0, then eulerian
    // odd count never is 1 for undirected graph

    return (odd == 2) ? 2 : 1;
}

igraph_bool_t check_if_bridge(igraph_t *g, int start, int v, long edge) {

    igraph_vector_t bridges;
    igraph_vector_init(&bridges, 0);
    igraph_bridges(g, &bridges);

    for (int i = 0; i < igraph_vector_size(&bridges); i++) {
        if (VECTOR(bridges)[i] == edge) {
            return 0;   
        }
    }

    return 1;
}

void print_euler_tool_implementation(int start, igraph_t *g) {

    int i = 0;
    while (true) {
        igraph_inclist_t il;
        igraph_inclist_init(g, &il, IGRAPH_ALL);
        igraph_vector_int_t *incedges;
        incedges = igraph_inclist_get(&il, start);
        long nc = igraph_vector_int_size(incedges);

        if (nc == 0) break;

        long edge = (long) VECTOR(*incedges)[i];
        igraph_integer_t v = IGRAPH_TO(g, edge) == start ? IGRAPH_FROM(g, edge) : IGRAPH_TO(g, edge);
        if (nc == 1 || check_if_bridge(g, start, v, edge)) {
            printf("%d - %d, ", start, v);
            igraph_es_t edge_to_delete;
            igraph_es_1(&edge_to_delete, edge);
            igraph_delete_edges(g, edge_to_delete);
            print_euler_tool_implementation(v, g);
        }
        i++;
    }

// the below implementation is inferior, as nc always changes

/*
    igraph_inclist_t il;
    igraph_inclist_init(g, &il, IGRAPH_ALL);
    igraph_vector_int_t *incedges;
    incedges = igraph_inclist_get(&il, start);
    long nc = igraph_vector_int_size(incedges);

    for (int i = 0; i < nc; i++) {
        long edge = (long) VECTOR(*incedges)[i];
        igraph_integer_t v = IGRAPH_TO(g, edge) == start ? IGRAPH_FROM(g, edge) : IGRAPH_TO(g, edge);
        if (nc == 1 || check_if_bridge(g, start, v, edge)) {
            printf("%d - %d, ", start, v);
            igraph_es_t edge_to_delete;
            igraph_es_1(&edge_to_delete, edge);
            igraph_delete_edges(g, edge_to_delete);
            print_euler_tool_implementation(v, g);
        }
    }
*/
}

void print_euler_tour(igraph_t *graph) {

    // default starting vertex is 0
    // when we have no odd vertices, we can start anywhere (usually 0)
    // when we have 2 odd vertices, we start at any one of them

    int res = is_eulerian(graph);

    if (res == 0) {
        IGRAPH_WARNING("Euler cycle not possible");
        return;
    }

    igraph_t copy = *graph; // making a copy of the graph
    // this is because this algorithm will delete edges 
    // therefore we do not want to modify our original graph

    int start = 0;
    //print_euler_tool_implementation(start, &copy);
    //printf("\n done\n");

    if (res == 1) {
        print_euler_tool_implementation(start, graph);
    } else {
        igraph_inclist_t il;
        igraph_inclist_init(graph, &il, IGRAPH_ALL);
        igraph_vector_int_t *incedges;
        for (int i = 0; i < igraph_vcount(graph); i++) {
            incedges = igraph_inclist_get(&il, i);
            if ( igraph_vector_int_size(incedges) % 2 == 1) {
                start = i;
                break;
            }
        }
        print_euler_tool_implementation(start, graph);
    }
    printf("\n done\n");
}

int main(void) {
    igraph_t graph;

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, -1);
    assert(is_eulerian(&graph) ==  2);
    print_euler_tour(&graph);

    printf("test 1 done\n");
    
    igraph_t graph2;
    igraph_small(&graph2, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, 5,2,
                2,6, 6,4, 4,8, 2,8, 2,7, 0,7, -1);
    assert(is_eulerian(&graph2) ==  1);
    print_euler_tour(&graph2);
    
    printf("test 2 done\n");

    igraph_t graph3;
    igraph_small(&graph3, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 2,4 , 3,5 , 4,5,
                4,6, 0,6, 6,7, 1,7, -1);
    assert(is_eulerian(&graph3) ==  0);

    printf("test 3 done\n");

    igraph_t graph4;
    igraph_small(&graph4, 0, IGRAPH_UNDIRECTED, 0,1 , 1,2, 2,3, 3,4 , 2,4 , 1,5,
                0,5 , -1);
    assert(is_eulerian(&graph4) ==  2);
    print_euler_tour(&graph4);

    printf("test 4 done\n");

    igraph_t graph5;
    igraph_small(&graph5, 0, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,4, 3,4, 1,3, 2,5, 4,5, 2,6, 1,6, 0,4, 6,5, -1);
    assert(is_eulerian(&graph5) == 2);
    print_euler_tour(&graph5);

    printf("test 5 done\n");

    printf("all tests done\n");


    return 0;
}