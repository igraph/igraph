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

#include "igraph_eulerian.h"
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

/* solution adapted from https://www.geeksforgeeks.org/eulerian-path-and-circuit/
The function returns one of the following values
0 -> if graph is not Eulerian
1 -> if graph has an euler circuit
2 -> if graph has an euler path */
int is_eulerian_undirected(igraph_t *graph, igraph_bool_t *has_path, igraph_bool_t *has_cycle) {

    igraph_bool_t res;
    igraph_integer_t odd, i;
    igraph_vector_t degree;
    
    igraph_vector_init(&degree, 0);

    odd = 0;

    IGRAPH_CHECK(igraph_is_connected(graph, &res, IGRAPH_WEAK));
    if (!res) {
        *has_path = 0;
        *has_cycle = 0;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_vector_init(&degree, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &degree);

    igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    for (i = 0; i < igraph_vector_size(&degree); i++) {
        if (((long int) VECTOR(degree)[i]) % 2 == 1) odd++;
    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    if (odd > 2) {
        *has_path = 0;
        *has_cycle = 0;
        return IGRAPH_SUCCESS;
    }

    /* if pdd count is 2, then semi-eulerian
    if odd count is 0, then eulerian
    odd count never is 1 for undirected graph
    */

    if (odd == 2) {
        *has_path = 1;
        *has_cycle = 0;
    } else {
        *has_path = 0;
        *has_cycle = 1;
    }

    return IGRAPH_SUCCESS;
}

int is_eulerian_directed(igraph_t *graph, igraph_bool_t *has_path, igraph_bool_t *has_cycle) {
    igraph_bool_t res_strong, res_weak;
    igraph_integer_t incoming_excess, outgoing_excess;
    igraph_integer_t i;
    igraph_vector_t out_degree, in_degree;

    incoming_excess = 0;
    outgoing_excess = 0;

    IGRAPH_CHECK(igraph_vector_init(&out_degree, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &out_degree);
    igraph_degree(graph, &out_degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);

    IGRAPH_CHECK(igraph_vector_init(&in_degree, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &in_degree);
    igraph_degree(graph, &in_degree, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);

    /* checking if incoming vertices == outgoing vertices */
    for (i = 0; i < igraph_vcount(graph); i++) {
        if (VECTOR(in_degree)[i] != VECTOR(out_degree)[i]) {
            if ((VECTOR(in_degree)[i] + 1 == VECTOR(out_degree)[i]) && (incoming_excess < 2 && outgoing_excess < 1)) {
                outgoing_excess++;
            } else if ((VECTOR(out_degree)[i] + 1 == VECTOR(in_degree)[i]) && (incoming_excess < 1 && outgoing_excess < 2)) {
                incoming_excess++;
            } else {
                *has_path = 0;
                *has_cycle = 0;
                igraph_vector_destroy(&in_degree);
                igraph_vector_destroy(&out_degree);
                return IGRAPH_SUCCESS;
            }
        }
    }

    igraph_vector_destroy(&in_degree);
    igraph_vector_destroy(&out_degree);

    IGRAPH_CHECK(igraph_is_connected(graph, &res_strong, IGRAPH_STRONG));
    IGRAPH_CHECK(igraph_is_connected(graph, &res_weak, IGRAPH_WEAK));

    if ((outgoing_excess == 0 && incoming_excess == 0) && (res_strong)) {
        *has_path = 0;
        *has_cycle = 1;
        return IGRAPH_SUCCESS;
    } else if ((outgoing_excess == 1 && incoming_excess == 1) && 
        (res_strong || res_weak)) {
        *has_path = 1;
        *has_cycle = 0;
        return IGRAPH_SUCCESS;
    }

    *has_path = 0;
    *has_cycle = 0;

    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/* flexible function */
int igraph_is_eulerian(igraph_t *graph, igraph_bool_t *has_path, igraph_bool_t *has_cycle) {
    if (igraph_is_directed(graph)) {
        is_eulerian_directed(graph, has_path, has_cycle);
    } else {
        is_eulerian_undirected(graph, has_path, has_cycle);
    } 
    return IGRAPH_SUCCESS;
}

igraph_bool_t check_if_bridge(igraph_t *g, igraph_integer_t start, igraph_integer_t v, long edge) {

    igraph_vector_t bridges;
    igraph_integer_t i;

    IGRAPH_CHECK(igraph_vector_init(&bridges, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &bridges);
    igraph_bridges(g, &bridges);

    for (i = 0; i < igraph_vector_size(&bridges); i++) {
        if (VECTOR(bridges)[i] == edge) {
            igraph_vector_destroy(&bridges);
            IGRAPH_FINALLY_CLEAN(1);
            return 0;   
        }
    }

    igraph_vector_destroy(&bridges);
    IGRAPH_FINALLY_CLEAN(1);

    return 1;
}

int print_euler_undirected_implementation(igraph_integer_t start, igraph_t *g, igraph_vector_t *path) {

    igraph_integer_t i, v;
    igraph_inclist_t il;
    igraph_vector_int_t *incedges;
    igraph_es_t edge_to_delete;
    long nc, edge;

    i = 0;

    while (1) {
        //igraph_inclist_init(g, &il, IGRAPH_ALL);
        IGRAPH_CHECK(igraph_inclist_init(g, &il, IGRAPH_ALL)); 
        IGRAPH_FINALLY(igraph_inclist_destroy, &il); 
        incedges = igraph_inclist_get(&il, start);
        nc = igraph_vector_int_size(incedges);

        if (nc == 0) break;

        edge = (long) VECTOR(*incedges)[i];
        v = IGRAPH_TO(g, edge) == start ? IGRAPH_FROM(g, edge) : IGRAPH_TO(g, edge);
        if (nc == 1 || check_if_bridge(g, start, v, edge)) {
            IGRAPH_CHECK(igraph_vector_push_back(path, v));
            IGRAPH_CHECK(igraph_es_1(&edge_to_delete, edge));
            IGRAPH_CHECK(igraph_delete_edges(g, edge_to_delete));
            print_euler_undirected_implementation(v, g, path);
        }
        i++;
    }


    igraph_inclist_destroy(&il);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;

}


/* solution adapted from https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/ */
int igraph_euler_path_undirected(igraph_t *graph, igraph_vector_t *path) {

    /* default starting vertex is 0
    when we have no odd vertices, we can start anywhere (usually 0)
    when we have 2 odd vertices, we start at any one of them
    igraph_vector_t path must be initialised prior to its us */

    igraph_integer_t start, i;
    igraph_t copy; 
    igraph_inclist_t incl;
    igraph_vector_int_t *incedges;
    igraph_bool_t has_path;
    igraph_bool_t cycle;

    has_path = 0;
    cycle = 0;

    igraph_is_eulerian(graph, &has_path, &cycle);

    if (!has_path && !cycle) {
        IGRAPH_WARNING("Euler cycle not possible");
        return IGRAPH_FAILURE;
    }

    IGRAPH_CHECK(igraph_copy(&copy, graph));
    IGRAPH_FINALLY(igraph_destroy, &copy);
    /* making a copy of the graph
    this is because this algorithm will delete edges 
    therefore we do not want to modify our original graph */

    start = 0;

    IGRAPH_CHECK(igraph_inclist_init(graph, &incl, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incl);

    if (cycle && !has_path) {
        IGRAPH_CHECK(igraph_vector_push_back(path, start));
        print_euler_undirected_implementation(start, &copy, path);
    } else {
        for (i = 0; i < igraph_vcount(graph); i++) {
            incedges = igraph_inclist_get(&incl, i);
            if ( igraph_vector_int_size(incedges) % 2 == 1) {
                start = i;
                break;
            }
        }
        IGRAPH_CHECK(igraph_vector_push_back(path, start));
        print_euler_undirected_implementation(start, &copy, path);
    }

    igraph_inclist_destroy(&incl);
    igraph_destroy(&copy);

    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;

}

/* solution adapted from https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/ */
int eulerian_path_directed_implementation(igraph_t *graph, igraph_integer_t *start_node, igraph_vector_t *outgoing_list, igraph_vector_t *res) {
    igraph_integer_t start = *start_node;
    igraph_integer_t curr = start;
    igraph_integer_t next;
    igraph_inclist_t il;
    igraph_es_t es;
    igraph_stack_t path;
    igraph_stack_t tracker;
    igraph_vector_int_t *incedges;
    long nc;
    long edge;

    IGRAPH_CHECK(igraph_stack_init(&path, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_stack_destroy, &path);

    IGRAPH_CHECK(igraph_stack_init(&tracker, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_stack_destroy, &tracker);

    igraph_stack_push(&tracker, start);

    while (!igraph_stack_empty(&tracker)) {
        
        if (VECTOR(*outgoing_list)[curr] != 0) {
            IGRAPH_CHECK(igraph_stack_push(&tracker, curr));
            
            IGRAPH_CHECK(igraph_inclist_init(graph, &il, IGRAPH_ALL));
            IGRAPH_FINALLY(igraph_inclist_destroy, &il);
            
            incedges = igraph_inclist_get(&il, curr);
            nc = igraph_vector_int_size(incedges);
            edge = (long) VECTOR(*incedges)[0];
            
            next = IGRAPH_TO(graph, edge);

            /* remove edge here */
            VECTOR(*outgoing_list)[curr]--;
            IGRAPH_CHECK(igraph_es_1(&es, edge));
            IGRAPH_CHECK(igraph_delete_edges(graph, es));

            igraph_inclist_destroy(&il);
            IGRAPH_FINALLY_CLEAN(1);

            curr = next;
        } else { /* back track to find remaining circuit */
            igraph_stack_push(&path, curr);
            curr = igraph_stack_pop(&tracker);
        }
    }

    while (!igraph_stack_empty(&path)) {
        IGRAPH_CHECK(igraph_vector_push_back(res, igraph_stack_pop(&path)));           
    }

    igraph_stack_destroy(&path);
    igraph_stack_destroy(&tracker);

    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


int igraph_eulerian_path_directed(igraph_t *graph, igraph_vector_t *res) {
    igraph_bool_t res_strong, res_weak;
    igraph_integer_t incoming_excess, outgoing_excess;
    igraph_integer_t start_node;
    igraph_t copy;
    igraph_integer_t i;
    igraph_vector_t out_degree, in_degree;

    IGRAPH_CHECK(igraph_copy(&copy, graph));
    IGRAPH_FINALLY(igraph_destroy, &copy);

    incoming_excess = 0;
    outgoing_excess = 0;

    /* determining the start node */
    /* also getting the outgoing list for each vector*/

    IGRAPH_CHECK(igraph_vector_init(&out_degree, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &out_degree);
    igraph_degree(graph, &out_degree, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);

    IGRAPH_CHECK(igraph_vector_init(&in_degree, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &in_degree);
    igraph_degree(graph, &in_degree, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);

    /* checking if incoming vertices == outgoing vertices */
    for (i = 0; i < igraph_vcount(graph); i++) {
        if (VECTOR(in_degree)[i] != VECTOR(out_degree)[i]) {
            if ((VECTOR(in_degree)[i] + 1 == VECTOR(out_degree)[i]) && (incoming_excess < 2 && outgoing_excess < 1)) {
                outgoing_excess++;
                start_node = i;
            } else if ((VECTOR(out_degree)[i] + 1 == VECTOR(in_degree)[i]) && (incoming_excess < 1 && outgoing_excess < 2)) {
                incoming_excess++;
            } 
        }
    }

    IGRAPH_CHECK(igraph_is_connected(graph, &res_strong, IGRAPH_STRONG));
    IGRAPH_CHECK(igraph_is_connected(graph, &res_weak, IGRAPH_WEAK));

    if ((outgoing_excess == 0 && incoming_excess == 0) && (res_strong)) {
        start_node = 0;
        eulerian_path_directed_implementation(&copy, &start_node, &out_degree, res);
    } else if ((outgoing_excess == 1 && incoming_excess == 1) && 
        (res_strong || res_weak)) {
        eulerian_path_directed_implementation(&copy, &start_node, &out_degree, res);
    }

    igraph_vector_destroy(&in_degree);
    igraph_vector_destroy(&out_degree);
    igraph_destroy(&copy);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


int igraph_eulerian_cycle(igraph_t *graph, igraph_vector_t *res) {
    double vec1, vec2;
    igraph_integer_t curr_edge;
    igraph_bool_t error;
    igraph_vector_t vector_res;
    igraph_bool_t cycle;
    igraph_bool_t has_path;

    has_path = 0;
    cycle = 0;

    igraph_is_eulerian(graph, &has_path, &cycle);

    if (!cycle) {
        IGRAPH_WARNING("Eulerian cycle not possible");
        return IGRAPH_FAILURE;
    }

    error = 0;
    IGRAPH_CHECK(igraph_vector_init(&vector_res, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector_res);

    if (igraph_is_directed(graph)) {
        igraph_eulerian_path_directed(graph, &vector_res);

        for (int i = 0; i < igraph_vector_size(&vector_res) - 1; i++) {
            vec1 = VECTOR(vector_res)[i];
            vec2 = VECTOR(vector_res)[i+1];

            igraph_get_eid(graph, &curr_edge, (int) vec1, (int) vec2, 1, error);
            igraph_vector_push_back(res, curr_edge);
        }
    } else {
        igraph_euler_path_undirected(graph, &vector_res);

        for (int i = 0; i < igraph_vector_size(&vector_res) - 1; i++) {
            vec1 = VECTOR(vector_res)[i];
            vec2 = VECTOR(vector_res)[i+1];

            igraph_get_eid(graph, &curr_edge, (int) vec1, (int) vec2, 0, error);
            igraph_vector_push_back(res, curr_edge);
        }
    } 
    igraph_vector_destroy(&vector_res);

    return IGRAPH_SUCCESS;
}

int igraph_eulerian_path(igraph_t *graph, igraph_vector_t *res) {
    double vec1, vec2;
    igraph_integer_t curr_edge;
    igraph_bool_t error;
    igraph_vector_t vector_res;
    igraph_bool_t cycle;
    igraph_bool_t has_path;

    has_path = 0;
    cycle = 0;

    igraph_is_eulerian(graph, &has_path, &cycle);

    if (!has_path && !cycle) {
        IGRAPH_WARNING("Eulerian path not possible");
        return IGRAPH_FAILURE;
    }

    error = 0;
    IGRAPH_CHECK(igraph_vector_init(&vector_res, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector_res);

    if (igraph_is_directed(graph)) {
        igraph_eulerian_path_directed(graph, &vector_res);

        for (int i = 0; i < igraph_vector_size(&vector_res) - 1; i++) {
            vec1 = VECTOR(vector_res)[i];
            vec2 = VECTOR(vector_res)[i+1];

            igraph_get_eid(graph, &curr_edge, (int) vec1, (int) vec2, 1, error);
            igraph_vector_push_back(res, curr_edge);
        }
    } else {
        igraph_euler_path_undirected(graph, &vector_res);

        for (int i = 0; i < igraph_vector_size(&vector_res) - 1; i++) {
            vec1 = VECTOR(vector_res)[i];
            vec2 = VECTOR(vector_res)[i+1];

            igraph_get_eid(graph, &curr_edge, (int) vec1, (int) vec2, 0, error);
            igraph_vector_push_back(res, curr_edge);
        }
    } 
    igraph_vector_destroy(&vector_res);

    return IGRAPH_SUCCESS;
}