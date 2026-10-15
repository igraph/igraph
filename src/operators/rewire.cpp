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

#include "core/exceptions.h"
#include "operators/rewire_internal.h"

#include <vector>
#include <unordered_map>

// TODO: fix current test errors
// after that replace with custom adjlist and see if it works
// ONLY HANDLE SIMPLE, UNDIRECTED right now

struct Edge {
    igraph_int_t id;
    igraph_int_t from;
    igraph_int_t to;

    bool operator==(const Edge& other) const {
        return from == other.from && to == other.to;
    }
};

struct AdjList {
    using NeighborMap = std::unordered_map<igraph_int_t, igraph_int_t>;
    igraph_int_t n_nodes;
    igraph_int_t n_edges;
    std::vector<NeighborMap> adjlist;

    explicit AdjList(const igraph_t *graph, const igraph_vector_int_t& edgelist)
        : n_nodes(igraph_vcount(graph)), adjlist(n_nodes)
    {
        // put edges into adjlist
        n_edges = igraph_ecount(graph);
        for (igraph_int_t i = 0; i < n_edges; i++) {
            igraph_int_t n1 = VECTOR(edgelist)[2 * i];
            igraph_int_t n2 = VECTOR(edgelist)[2 * i + 1];

            // undirected implementation
            // TODO: if we have loops later this needs to change (double counts self loops) 
            add_edge(n1, n2);
        }
    }

    void add_edge(igraph_int_t n1, igraph_int_t n2) {
        adjlist[n1][n2]++;
        adjlist[n2][n1]++;
    }

    void remove_edge(igraph_int_t n1, igraph_int_t n2) {
        adjlist[n1][n2]--;
        adjlist[n2][n1]--;

        if (adjlist[n1][n2] <= 0) { // should be symmetrical
            adjlist[n1].erase(n2);  // should be O(1) average
            adjlist[n2].erase(n1);
        }
    }

    // function: has_edge()
    bool has_edge(igraph_int_t n1, igraph_int_t n2) {
        // TODO: maybe add a check that what's found does not have edge <= 0
        // for robustness
        return adjlist[n1].find(n2) != adjlist[n1].end();
    }

    // function: replace_edge()
    // undirected only
    void replace_edge(igraph_int_t from, igraph_int_t old_to, igraph_int_t new_to) {
        remove_edge(from, old_to);
        add_edge(from, new_to);
    }
    // make sure to delete keys that go to 0, has_edge() depends on it
};

// function to pick random edges
Edge get_random_edge(const igraph_vector_int_t& edgelist, igraph_int_t n_edges) {
    // get random edge
    igraph_int_t eid = RNG_INTEGER(0, n_edges - 1); // high inclusive

    // get endpoints of edge
    igraph_int_t from = VECTOR(edgelist)[2 * eid];
    igraph_int_t to = VECTOR(edgelist)[2 * eid + 1];

    return Edge{eid, from, to};
}

// TODO: update when expanding beyond simple graphs
// does not destroy, but also does not create, multiedges
// loops later
igraph_bool_t can_rewire(Edge e1, Edge e2, AdjList& adjlist) {
    if (e1.from == e2.to || e1.to == e2.from) return false; // loop
    // TODO: discussed rewiring anyway, but then adjlist_replace_edge throws error because
    // edge already exists
    else if (e1.from == e2.from || e1.to == e2.to) return false; // no effect
    // TODO: specifically stop rewiring if we picked the EXACT same edge?
    
    // check if results multiedge
    if (adjlist.has_edge(e1.from, e2.to)) return false;
    if (adjlist.has_edge(e2.from, e1.to)) return false;

    return true;
}

// attempts to perform one rewiring
igraph_error_t rewire_once(igraph_vector_int_t& edgelist, AdjList& adjlist, 
                          igraph_int_t n_edges) {
    // pick edges from EDGELIST
    Edge e1 = get_random_edge(edgelist, n_edges);
    Edge e2 = get_random_edge(edgelist, n_edges);

    // with .5 probability, flip edge 2 (SIMPLE GRAPHS)
    if (RNG_BOOL()) {
        igraph_int_t temp = e2.from; e2.from = e2.to; e2.to = temp;
    }
    
    if (can_rewire(e1, e2, adjlist)) {
        // update edgelist
        VECTOR(edgelist)[e1.id * 2 + 1] = e2.to;
        VECTOR(edgelist)[e2.id * 2] = e2.from; // edge 2 is potentially flipped
        VECTOR(edgelist)[e2.id * 2 + 1] = e1.to;

        // update adjlist
        // TODO: same thing with directed param
        adjlist.replace_edge(e1.from, e1.to, e2.to);
        adjlist.replace_edge(e2.from, e2.to, e1.to);

    } else {
        // TODO: track failures/no effects too
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_rewire_2(igraph_t *graph, igraph_int_t n, 
                               igraph_edge_type_sw_t allowed_edge_types) {

    const igraph_int_t n_edges = igraph_ecount(graph);

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    // set up edgelist
    igraph_vector_int_t edgelist;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edgelist, n_edges * 2);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edgelist, /*bycol=*/ false));

    // set up adjlist
    AdjList adjlist(graph, edgelist);

    for (igraph_int_t i = 0; i < n; i++) {
        // TODO: keep track of % successful rewirings
        IGRAPH_CHECK(rewire_once(edgelist, adjlist, n_edges));
    }

    // write changes to graph
    IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID)));
    IGRAPH_CHECK(igraph_add_edges(graph, &edgelist, 0));

    // free memory and stuff
    igraph_vector_int_destroy(&edgelist);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}
