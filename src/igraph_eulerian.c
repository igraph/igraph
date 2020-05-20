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
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_components.h"
#include "igraph_types_internal.h"

/* solution adapted from https://www.geeksforgeeks.org/eulerian-path-and-circuit/
The function returns one of the following values
has_path is set to 1 if a path exists, 0 otherwise
has_cycle is set to 1 if a cycle exists, 0 otherwise
*/
int is_eulerian_undirected(igraph_t *graph, igraph_bool_t *has_path, igraph_bool_t *has_cycle) {

    igraph_bool_t res;
    igraph_integer_t odd, i;
    igraph_vector_t degree;

    IGRAPH_CHECK(igraph_is_connected(graph, &res, IGRAPH_WEAK));
    if (!res) {
        *has_path = 0;
        *has_cycle = 0;
        return IGRAPH_SUCCESS;
    }

    odd = 0;

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
                IGRAPH_FINALLY_CLEAN(2);

                return IGRAPH_SUCCESS;
            }
        }
    }

    igraph_vector_destroy(&in_degree);
    igraph_vector_destroy(&out_degree);

    IGRAPH_FINALLY_CLEAN(2);

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

    return IGRAPH_SUCCESS;
}

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

    igraph_vector_init(&bridges, 0);
    igraph_bridges(g, &bridges);

    for (i = 0; i < igraph_vector_size(&bridges); i++) {
        if (VECTOR(bridges)[i] == edge) {
            igraph_vector_destroy(&bridges);
            return 0;   
        }
    }

    igraph_vector_destroy(&bridges);

    return 1;
}

int euler_undirected_implementation(igraph_integer_t start, igraph_t *g, igraph_vector_t *path, igraph_vector_bool_t *visited_list) {

    igraph_integer_t i, v, adj;
    igraph_inclist_t il;
    igraph_vector_int_t *incedges;
    long nc, edge;
    int j;

    i = 0;

    while (1) {
        IGRAPH_CHECK(igraph_inclist_init(g, &il, IGRAPH_ALL)); 
        IGRAPH_FINALLY(igraph_inclist_destroy, &il); 
        incedges = igraph_inclist_get(&il, start);
        nc = igraph_vector_int_size(incedges);
        adj = 0;

        for (j = 0; j < nc; j++) {
            edge = (long) VECTOR(*incedges)[j];
            if (!VECTOR(*visited_list)[edge]) adj++;
        }

        if (adj == 0) break;

        edge = (long) VECTOR(*incedges)[i];
        if (!VECTOR(*visited_list)[edge]) {
            v = IGRAPH_TO(g, edge) == start ? IGRAPH_FROM(g, edge) : IGRAPH_TO(g, edge);
            if (adj == 1 || check_if_bridge(g, start, v, edge)) {
                VECTOR(*visited_list)[edge] = 1;
                IGRAPH_CHECK(igraph_vector_push_back(path, edge));               
                euler_undirected_implementation(v, g, path, visited_list);
            }
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
    when we have 2 odd vertices, we start at any one of them */

    igraph_integer_t start, i, e_count;
    igraph_inclist_t incl;
    igraph_vector_int_t *incedges;
    igraph_bool_t has_path;
    igraph_bool_t cycle;
    igraph_vector_bool_t visited_list;

    has_path = 0;
    cycle = 0;

    igraph_is_eulerian(graph, &has_path, &cycle);

    e_count = igraph_ecount(graph);
    igraph_vector_bool_init(&visited_list, e_count);

    for (i = 0; i < e_count; i++) {
        VECTOR(visited_list)[i] = 0;
    }

    start = 0;

    IGRAPH_CHECK(igraph_inclist_init(graph, &incl, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_inclist_destroy, &incl);

    if (cycle && !has_path) {
        euler_undirected_implementation(start, graph, path, &visited_list);
    } else {
        for (i = 0; i < igraph_vcount(graph); i++) {
            incedges = igraph_inclist_get(&incl, i);
            if ( igraph_vector_int_size(incedges) % 2 == 1) {
                start = i;
                break;
            }
        }
        euler_undirected_implementation(start, graph, path, &visited_list);
    }

    igraph_inclist_destroy(&incl);

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;

}

/* solution adapted from https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/ */
int eulerian_path_directed_implementation(igraph_t *graph, igraph_integer_t *start_node, igraph_vector_t *outgoing_list, igraph_vector_t *res) {
    igraph_integer_t start = *start_node;
    igraph_integer_t curr = start;
    igraph_integer_t next, e_count, curr_e;
    igraph_inclist_t il;
    igraph_stack_t path, tracker, edge_tracker, edge_path;
    igraph_vector_int_t *incedges;
    igraph_vector_bool_t visited_list;
    long nc, edge;
    int i, j;

    IGRAPH_CHECK(igraph_stack_init(&path, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_stack_destroy, &path);

    IGRAPH_CHECK(igraph_stack_init(&tracker, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_stack_destroy, &tracker);

    IGRAPH_CHECK(igraph_stack_init(&edge_path, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_stack_destroy, &edge_path);

    IGRAPH_CHECK(igraph_stack_init(&edge_tracker, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_stack_destroy, &edge_tracker);

    e_count = igraph_ecount(graph);
    igraph_vector_bool_init(&visited_list, e_count);

    for (i = 0; i < e_count; i++) {
        VECTOR(visited_list)[i] = 0;
    }

    igraph_stack_push(&tracker, start);

    while (!igraph_stack_empty(&tracker)) {
        
        if (VECTOR(*outgoing_list)[curr] != 0) {
            IGRAPH_CHECK(igraph_stack_push(&tracker, curr));
            
            IGRAPH_CHECK(igraph_inclist_init(graph, &il, IGRAPH_OUT));
            IGRAPH_FINALLY(igraph_inclist_destroy, &il);
            
            incedges = igraph_inclist_get(&il, curr);
            nc = igraph_vector_int_size(incedges);

            for (j = 0; j < nc; j++) {
                edge = (long) VECTOR(*incedges)[j];
                if (!VECTOR(visited_list)[edge]) {
                    break;
                }
            }
            
            next = IGRAPH_TO(graph, edge);

            IGRAPH_CHECK(igraph_stack_push(&edge_tracker, edge));

            /* remove edge here */
            VECTOR(*outgoing_list)[curr]--;
            VECTOR(visited_list)[edge] = 1;
            
            igraph_inclist_destroy(&il);
            IGRAPH_FINALLY_CLEAN(1);

            curr = next;
        } else { /* back track to find remaining circuit */
            igraph_stack_push(&path, curr);
            curr = igraph_stack_pop(&tracker);
            if (!igraph_stack_empty(&edge_tracker)) {
                curr_e = igraph_stack_pop(&edge_tracker);
                igraph_stack_push(&edge_path, curr_e);
            }
        }
    }

    while (!igraph_stack_empty(&edge_path)) {
        IGRAPH_CHECK(igraph_vector_push_back(res, igraph_stack_pop(&edge_path)));           
    }

    igraph_stack_destroy(&path);
    igraph_stack_destroy(&tracker);
    igraph_stack_destroy(&edge_path);
    igraph_stack_destroy(&edge_tracker);

    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


int igraph_eulerian_path_directed(igraph_t *graph, igraph_vector_t *res) {
    igraph_bool_t res_strong, res_weak;
    igraph_integer_t incoming_excess, outgoing_excess;
    igraph_integer_t start_node;
    igraph_integer_t i;
    igraph_vector_t out_degree, in_degree;

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
        eulerian_path_directed_implementation(graph, &start_node, &out_degree, res);
    } else if ((outgoing_excess == 1 && incoming_excess == 1) && 
        (res_strong || res_weak)) {
        eulerian_path_directed_implementation(graph, &start_node, &out_degree, res);
    }

    igraph_vector_destroy(&in_degree);
    igraph_vector_destroy(&out_degree);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


int igraph_eulerian_cycle(igraph_t *graph, igraph_vector_t *res) {
    igraph_bool_t cycle;
    igraph_bool_t has_path;

    has_path = 0;
    cycle = 0;

    igraph_is_eulerian(graph, &has_path, &cycle);

    if (!cycle) {
        return IGRAPH_EINVAL;
    }

    if (igraph_is_directed(graph)) {
        igraph_eulerian_path_directed(graph, res);
    } else {
        igraph_euler_path_undirected(graph, res);
    }

    return IGRAPH_SUCCESS;
}

int igraph_eulerian_path(igraph_t *graph, igraph_vector_t *res) {
    igraph_bool_t cycle;
    igraph_bool_t has_path;

    has_path = 0;
    cycle = 0;

    igraph_is_eulerian(graph, &has_path, &cycle);

    if (!has_path && !cycle) {
        return IGRAPH_EINVAL;
    }

    if (igraph_is_directed(graph)) {
        igraph_eulerian_path_directed(graph, res);
    } else {
        igraph_euler_path_undirected(graph, res);
    }

    return IGRAPH_SUCCESS;
}