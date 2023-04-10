/*
 *
 * gengraph - generation of random simple connected graphs with prescribed
 *            degree sequence
 *
 * Copyright (C) 2006  Fabien Viger
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GRAPH_MOLLOY_HASH_H
#define GRAPH_MOLLOY_HASH_H

#include "gengraph_definitions.h"
#include "gengraph_hash.h"
#include "gengraph_degree_sequence.h"

#include "igraph_datatype.h"

#include <string.h>
#include <assert.h>
// This class handles graphs with a constant degree sequence.

#define FINAL_HEURISTICS        0
#define GKAN_HEURISTICS         1
#define FAB_HEURISTICS          2
#define OPTIMAL_HEURISTICS      3
#define BRUTE_FORCE_HEURISTICS  4

namespace gengraph {

//****************************
//  class graph_molloy_hash
//****************************

class graph_molloy_hash {

private:
    // Number of vertices
    igraph_integer_t n;
    //Number of arcs ( = #edges * 2 )
    igraph_integer_t a;
    //Total size of links[]
    igraph_integer_t size;
    // The degree sequence of the graph
    igraph_integer_t *deg;
    // The array containing all links
    igraph_integer_t *links;
    // The array containing pointers to adjacency list of every vertices
    igraph_integer_t **neigh;
    // Counts total size
    void compute_size();
    // Build neigh with deg and links
    void compute_neigh();
    // Allocate memory according to degree_sequence (for constructor use only!!)
    igraph_integer_t alloc(degree_sequence &);
    // Add edge (u,v). Return FALSE if vertex a is already full.
    // WARNING : only to be used by havelhakimi(), restore() or constructors
    inline bool add_edge(igraph_integer_t u, igraph_integer_t v, igraph_integer_t *realdeg) {
        igraph_integer_t deg_u = realdeg[u];
        if (deg_u == deg[u]) {
            return false;
        }
        // Check that edge was not already inserted
        assert(fast_search(neigh[u], (u == n - 1 ? links + size : neigh[u + 1]) - neigh[u], v) == NULL);
        assert(fast_search(neigh[v], (v == n - 1 ? links + size : neigh[v + 1]) - neigh[v], u) == NULL);
        assert(deg[u] < deg_u);
        igraph_integer_t deg_v = realdeg[v];
        if (IS_HASH(deg_u)) {
            *H_add(neigh[u], HASH_EXPAND(deg_u), v) = v;
        } else {
            neigh[u][deg[u]] = v;
        }
        if (IS_HASH(deg_v)) {
            *H_add(neigh[v], HASH_EXPAND(deg_v), u) = u;
        } else {
            neigh[v][deg[v]] = u;
        }
        deg[u]++;
        deg[v]++;
        // Check that edge was actually inserted
        assert(fast_search(neigh[u], int((u == n - 1 ? links + size : neigh[u + 1]) - neigh[u]), v) != NULL);
        assert(fast_search(neigh[v], int((v == n - 1 ? links + size : neigh[v + 1]) - neigh[v]), u) != NULL);
        return true;
    }
    // Swap edges
    inline void swap_edges(igraph_integer_t from1, igraph_integer_t to1, igraph_integer_t from2, igraph_integer_t to2) {
        H_rpl(neigh[from1], deg[from1], to1, to2);
        H_rpl(neigh[from2], deg[from2], to2, to1);
        H_rpl(neigh[to1], deg[to1], from1, from2);
        H_rpl(neigh[to2], deg[to2], from2, from1);
    }
    // Backup graph [sizeof(igraph_integer_t) bytes per edge]
    igraph_integer_t* backup();
    // Test if vertex is in an isolated component of size<K
    bool isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited);
    // Pick random edge, and gives a corresponding vertex
    inline igraph_integer_t pick_random_vertex() {
        igraph_integer_t v;
        do {
            v = links[my_random() % size];
        } while (v == HASH_NONE);
        return v;
    }
    // Pick random neighbour
    inline igraph_integer_t* random_neighbour(igraph_integer_t v) {
        return H_random(neigh[v], deg[v]);
    }
    // Depth-first search.
    igraph_integer_t depth_search(bool *visited, igraph_integer_t *buff, igraph_integer_t v0 = 0);
    // Returns complexity of isolation test
    igraph_integer_t effective_isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited);
    // Depth-Exploration. Returns number of steps done. Stops when encounter vertex of degree > dmax.
    void depth_isolated(igraph_integer_t v, igraph_integer_t &calls, igraph_integer_t &left_to_explore, igraph_integer_t dmax, igraph_integer_t * &Kbuff, bool *visited);


public:
    //degree of v
    inline igraph_integer_t degree(igraph_integer_t v) {
        return deg[v];
    };
    // For debug purposes : verify validity of the graph (symetry, simplicity)
    //bool verify();
    // Destroy deg[], neigh[] and links[]
    ~graph_molloy_hash();
    // Allocate memory for the graph. Create deg and links. No edge is created.
    graph_molloy_hash(degree_sequence &);
    // Create graph from hard copy
    graph_molloy_hash(igraph_integer_t *);
    // Create hard copy of graph
    igraph_integer_t *hard_copy();
    // Restore from backup
    void restore(igraph_integer_t* back);
    //Clear hash tables
    void init();
    // nb arcs
    inline igraph_integer_t nbarcs() {
        return a;
    };
    // nb vertices
    inline igraph_integer_t nbvertices() {
        return n;
    };
    // print graph in SUCC_LIST mode, in stdout
    /* void print(FILE *f = stdout); */
    igraph_error_t print(igraph_t *graph);
    // Test if graph is connected
    bool is_connected();
    // is edge ?
    inline bool is_edge(igraph_integer_t u, igraph_integer_t v) {
        assert(H_is(neigh[u], deg[u], v) == (fast_search(neigh[u], HASH_SIZE(deg[u]), v) != NULL));
        assert(H_is(neigh[v], deg[v], u) == (fast_search(neigh[v], HASH_SIZE(deg[v]), u) != NULL));
        assert(H_is(neigh[u], deg[u], v) == H_is(neigh[v], deg[v], u));
        if (deg[u] < deg[v]) {
            return H_is(neigh[u], deg[u], v);
        } else {
            return H_is(neigh[v], deg[v], u);
        }
    }
    // Random edge swap ATTEMPT. Return 1 if attempt was a succes, 0 otherwise
    int random_edge_swap(igraph_integer_t K = 0, igraph_integer_t *Kbuff = NULL, bool *visited = NULL);
    // Connected Shuffle
    igraph_integer_t shuffle(igraph_integer_t, igraph_integer_t, int type);
    // Optimal window for the gkantsidis heuristics
    igraph_integer_t optimal_window();
    // Average unitary cost per post-validated edge swap, for some window
    double average_cost(igraph_integer_t T, igraph_integer_t *back, double min_cost);
    // Get caracteristic K
    double eval_K(int quality = 100);
    // Get effective K
    double effective_K(igraph_integer_t K, int quality = 10000);
    // Try to shuffle T times. Return true if at the end, the graph was still connected.
    bool try_shuffle(igraph_integer_t T, igraph_integer_t K, igraph_integer_t *back = NULL);
};

} // namespace gengraph

#endif //GRAPH_MOLLOY_HASH_H
