/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_operators.h"

#include "igraph_random.h"
#include "igraph_error.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"

#include "operators/rewire_internal.h"

struct Edge {
    igraph_integer_t id;
    igraph_integer_t from;
    igraph_integer_t to;

    bool operator==(const Edge& other) const {
        return from == other.from && to == other.to;
    }
};

// function to pick random edges
Edge get_random_edge(const igraph_vector_int_t& edgelist, igraph_integer_t n_edges) {
    // get random edge
    igraph_integer_t eid = RNG_INTEGER(0, n_edges - 1); // high inclusive

    // get endpoints of edge
    igraph_integer_t from = VECTOR(edgelist)[2 * eid];
    igraph_integer_t to = VECTOR(edgelist)[2 * eid + 1];

    return Edge{eid, from, to};
}

// TODO: update when expanding beyond simple graphs
// does not destroy, but also does not create, multiedges
// loops later
igraph_bool_t can_rewire(Edge e1, Edge e2, igraph_adjlist_t& adjlist) {
    if (e1.from == e2.to || e1.to == e2.from) return false; // loop
    // TODO: discussed rewiring anyway, but then adjlist_replace_edge throws error because
    // edge already exists
    else if (e1.from == e2.from || e1.to == e2.to) return false; // no effect
    // TODO: specifically stop rewiring if we picked the EXACT same edge?
    
    // check if resulting multiedge
    // TODO: add graph type later
    if (igraph_adjlist_has_edge(&adjlist, e1.from, e2.to, false)) return false;
    if (igraph_adjlist_has_edge(&adjlist, e2.from, e1.to, false)) return false;
    return true;
}

// attempts to perform one rewiring
igraph_error_t rewire_once(igraph_vector_int_t& edgelist, igraph_adjlist_t& adjlist, 
                          igraph_integer_t n_edges) {
    // pick edges from EDGELIST
    Edge e1 = get_random_edge(edgelist, n_edges);
    Edge e2 = get_random_edge(edgelist, n_edges);

    // with .5 probability, flip edge 2 (SIMPLE GRAPHS)
    if (RNG_BOOL()) {
        igraph_integer_t temp = e2.from; e2.from = e2.to; e2.to = temp;
    }
    
    if (can_rewire(e1, e2, adjlist)) {
        // update edgelist
        VECTOR(edgelist)[e1.id * 2 + 1] = e2.to;
        VECTOR(edgelist)[e2.id * 2] = e2.from; // edge 2 is potentially flipped
        VECTOR(edgelist)[e2.id * 2 + 1] = e1.to;

        // update adjlist
        // TODO: same thing with directed param
        IGRAPH_CHECK(igraph_adjlist_replace_edge(&adjlist, e1.from, e1.to, e2.to, false));
        IGRAPH_CHECK(igraph_adjlist_replace_edge(&adjlist, e2.from, e2.to, e1.to, false));
    } else {
        // TODO: track failures/no effects too
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_rewire_2(igraph_t *graph, igraph_integer_t n, 
                               igraph_edge_type_sw_t allowed_edge_types) {
    const igraph_integer_t n_edges = igraph_ecount(graph);

    // set up edgelist
    igraph_vector_int_t edgelist;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist, n_edges * 2);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, /*bycol=*/ false));

    // set up adjlist
    igraph_adjlist_t adjlist;
    // TODO: for now, simple graphs only. loops and multi parameters reconsider once ready
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    for (igraph_integer_t i = 0; i < n; i++) {
        // TODO: keep track of % successful rewirings
        IGRAPH_CHECK(rewire_once(edgelist, adjlist, n_edges));
    }

    // write changes to graph
    IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID)));
    IGRAPH_CHECK(igraph_add_edges(graph, &edgelist, 0));

    // free memory and stuff
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&edgelist);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}