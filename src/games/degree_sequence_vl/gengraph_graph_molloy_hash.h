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
    int n;
    //Number of arcs ( = #edges * 2 )
    int a;
    //Total size of links[]
    int size;
    // The degree sequence of the graph
    int *deg;
    // The array containing all links
    int *links;
    // The array containing pointers to adjacency list of every vertices
    int **neigh;
    // Counts total size
    void compute_size();
    // Build neigh with deg and links
    void compute_neigh();
    // Allocate memory according to degree_sequence (for constructor use only!!)
    int alloc(degree_sequence &);
    // Add edge (u,v). Return FALSE if vertex a is already full.
    // WARNING : only to be used by havelhakimi(), restore() or constructors
    inline bool add_edge(int u, int v, int *realdeg) {
        int deg_u = realdeg[u];
        if (deg_u == deg[u]) {
            return false;
        }
        // Check that edge was not already inserted
        assert(fast_search(neigh[u], int((u == n - 1 ? links + size : neigh[u + 1]) - neigh[u]), v) == NULL);
        assert(fast_search(neigh[v], int((v == n - 1 ? links + size : neigh[v + 1]) - neigh[v]), u) == NULL);
        assert(deg[u] < deg_u);
        int deg_v = realdeg[v];
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
    inline void swap_edges(int from1, int to1, int from2, int to2) {
        H_rpl(neigh[from1], deg[from1], to1, to2);
        H_rpl(neigh[from2], deg[from2], to2, to1);
        H_rpl(neigh[to1], deg[to1], from1, from2);
        H_rpl(neigh[to2], deg[to2], from2, from1);
    }
    // Backup graph [sizeof(int) bytes per edge]
    int* backup();
    // Test if vertex is in an isolated component of size<K
    bool isolated(int v, int K, int *Kbuff, bool *visited);
    // Pick random edge, and gives a corresponding vertex
    inline int pick_random_vertex() {
        int v;
        do {
            v = links[my_random() % size];
        } while (v == HASH_NONE);
        return v;
    }
    // Pick random neighbour
    inline int* random_neighbour(const int v) {
        return H_random(neigh[v], deg[v]);
    }
    // Depth-first search.
    int depth_search(bool *visited, int *buff, int v0 = 0);
    // Returns complexity of isolation test
    long effective_isolated(int v, int K, int *Kbuff, bool *visited);
    // Depth-Exploration. Returns number of steps done. Stops when encounter vertex of degree > dmax.
    void depth_isolated(int v, long &calls, int &left_to_explore, int dmax, int * &Kbuff, bool *visited);


public:
    //degree of v
    inline int degree(const int v) {
        return deg[v];
    };
    // For debug purposes : verify validity of the graph (symetry, simplicity)
    bool verify();
    // Destroy deg[], neigh[] and links[]
    ~graph_molloy_hash();
    // Allocate memory for the graph. Create deg and links. No edge is created.
    graph_molloy_hash(degree_sequence &);
    // Create graph from hard copy
    graph_molloy_hash(int *);
    // Create hard copy of graph
    int *hard_copy();
    // Restore from backup
    void restore(int* back);
    //Clear hash tables
    void init();
    // nb arcs
    inline int nbarcs() {
        return a;
    };
    // nb vertices
    inline int nbvertices() {
        return n;
    };
    // print graph in SUCC_LIST mode, in stdout
    void print(FILE *f = stdout);
    int print(igraph_t *graph);
    // Test if graph is connected
    bool is_connected();
    // is edge ?
    inline bool is_edge(int u, int v) {
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
    int random_edge_swap(int K = 0, int *Kbuff = NULL, bool *visited = NULL);
    // Connected Shuffle
    unsigned long shuffle(unsigned long, unsigned long, int type);
    // Optimal window for the gkantsidis heuristics
    int optimal_window();
    // Average unitary cost per post-validated edge swap, for some window
    double average_cost(int T, int *back, double min_cost);
    // Get caracteristic K
    double eval_K(int quality = 100);
    // Get effective K
    double effective_K(int K, int quality = 10000);
    // Try to shuffle T times. Return true if at the end, the graph was still connected.
    bool try_shuffle(int T, int K, int *back = NULL);


    /*_____________________________________________________________________________
      Not to use anymore : use graph_molloy_opt class instead

    private:
      // breadth-first search. Store the distance (modulo 3)  in dist[]. Returns eplorated component size.
      int width_search(unsigned char *dist, int *buff, int v0=0);

    public:
      // Create graph
      graph_molloy_hash(FILE *f);
      // Bind the graph avoiding multiple edges or self-edges (return false if fail)
      bool havelhakimi();
      // Get the graph connected  (return false if fail)
      bool make_connected();
      // "Fab" Shuffle (Optimized heuristic of Gkantsidis algo.)
      long long fab_connected_shuffle(long long);
      // Naive Shuffle
      long long slow_connected_shuffle(long long);
      // Maximum degree
      int max_degree();
      // compute vertex betweenness : for each vertex, a unique random shortest path is chosen.
      // this choice is consistent (if shortest path from a to c goes through b and then d,
      // then shortest path from a to d goes through b). If(trivial path), also count all the
      // shortest paths where vertex is an extremity
      int *vertex_betweenness_rsp(bool trivial_path);
      // same, but when multiple shortest path are possible, average the weights.
      double *vertex_betweenness_asp(bool trivial_path);
    //___________________________________________________________________________________
    */

};

} // namespace gengraph

#endif //GRAPH_MOLLOY_HASH_H

