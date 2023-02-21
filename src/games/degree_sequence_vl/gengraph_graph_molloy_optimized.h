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
#ifndef GRAPH_MOLLOY_OPT_H
#define GRAPH_MOLLOY_OPT_H

#include "gengraph_definitions.h"
#include "gengraph_degree_sequence.h"

#include <assert.h>
#include "gengraph_random.h"

namespace gengraph {

// This class handles graphs with a constant degree sequence.

class graph_molloy_opt {

private:
    // Random generator
    KW_RNG::RNG rng;
    // Number of vertices
    igraph_integer_t n;
    //Number of arcs ( = #edges * 2 )
    igraph_integer_t a;
    // The degree sequence of the graph
    igraph_integer_t *deg;
    // The array containing all links
    igraph_integer_t *links;
    // The array containing pointers to adjacency list of every vertices
    igraph_integer_t **neigh;
    // Allocate memory according to degree_sequence (for constructor use only!!)
    void alloc(degree_sequence &);
    // Compute #edges
    inline void refresh_nbarcs() {
        a = 0;
        for (igraph_integer_t* d = deg + n; d != deg; ) {
            a += *(--d);
        }
    }
    // Build neigh with deg and links
    void compute_neigh();
    // Swap edges. The swap MUST be valid !!!
    inline void swap_edges(igraph_integer_t from1, igraph_integer_t to1, igraph_integer_t from2, igraph_integer_t to2) {
        fast_rpl(neigh[from1], to1, to2);
        fast_rpl(neigh[from2], to2, to1);
        fast_rpl(neigh[to1], from1, from2);
        fast_rpl(neigh[to2], from2, from1);
    }

    // Swap edges only if they are simple. return false if unsuccessful.
    bool swap_edges_simple(igraph_integer_t, igraph_integer_t, igraph_integer_t, igraph_integer_t);
    // Test if vertex is in an isolated component of size<K
    bool isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited);
    // Pick random edge, and gives a corresponding vertex
    inline igraph_integer_t pick_random_vertex() {
        return links[my_random() % a];
    };
    // Pick random neighbour
    inline igraph_integer_t* random_neighbour(const igraph_integer_t v) {
        return neigh[v] + (my_random() % deg[v]);
    };
    // Returns complexity of isolation test
    igraph_integer_t effective_isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited);
    // Depth-Exploration. Returns number of steps done. Stops when encounter vertex of degree > dmax.
    void depth_isolated(igraph_integer_t v, igraph_integer_t &calls, igraph_integer_t &left_to_explore, igraph_integer_t dmax, igraph_integer_t * &Kbuff, bool *visited);
    // breadth-first search. Store the distance (modulo 3)  in dist[]. Returns eplorated component size.
    igraph_integer_t width_search(unsigned char *dist, igraph_integer_t *buff, igraph_integer_t v0 = 0, igraph_integer_t toclear = -1);
    // depth-first search.
    igraph_integer_t depth_search(bool *visited, igraph_integer_t *buff, igraph_integer_t v0 = 0);
    // breadth-first search that count the number of shortest paths going from src to each vertex
    igraph_integer_t breadth_path_search(igraph_integer_t src, igraph_integer_t *buff, double *paths, unsigned char *dist);
    // Return component indexes where vertices belong to, starting from 0,
    // sorted by size (biggest component has index 0)
    igraph_integer_t *components(igraph_integer_t *comp = NULL);

public:
    // neigh[]
    inline igraph_integer_t** neighbors() {
        return neigh;
    };
    // deg[]
    inline igraph_integer_t* degrees() {
        return deg;
    };
    //adjacency list of v
    inline igraph_integer_t* operator[](const igraph_integer_t v) {
        return neigh[v];
    };
    //degree of v
    inline igraph_integer_t degree(const igraph_integer_t v) {
        return deg[v];
    };
    // Detach deg[] and neigh[]
    void detach();
    // Destroy deg and links
    ~graph_molloy_opt();
    // Create graph from file (stdin not supported unless rewind() possible)
    //graph_molloy_opt(FILE *f);
    // Allocate memory for the graph. Create deg and links. No edge is created.
    graph_molloy_opt(degree_sequence &);
    // Create graph from hard copy
    graph_molloy_opt(igraph_integer_t *);
    // Create hard copy of graph
    igraph_integer_t *hard_copy();
    // Remove unused edges, updates neigh[], recreate links[]
    void clean();
    // nb arcs
    inline igraph_integer_t nbarcs() {
        return a;
    };
    // last degree
    inline igraph_integer_t last_degree() {
        return deg[n - 1];
    };
    // nb vertices
    inline igraph_integer_t nbvertices() {
        return n;
    };
    // nb vertices having degree > 0
    inline igraph_integer_t nbvertices_real() {
        igraph_integer_t s = 0;
        for (igraph_integer_t *d = deg + n; d-- != deg; ) {
            if (*d) {
                s++;
            }
        }
        return s;
    };
    // return list of vertices with degree > 0. Compute #vertices, if not given.
    igraph_integer_t *vertices_real(igraph_integer_t &nb_v);
    // print graph in SUCC_LIST mode, in stdout
    void print(FILE *f = stdout, bool NOZERO = true);
    // Bind the graph avoiding multiple edges or self-edges (return false if fail)
    bool havelhakimi();
    // Get the graph connected  (return false if fail)
    bool make_connected();
    // Test if graph is connected
    bool is_connected();
    // Maximum degree
    igraph_integer_t max_degree();
    // is edge ?
    inline bool is_edge(const igraph_integer_t u, const igraph_integer_t v) {
        if (deg[v] < deg[u]) {
            return (fast_search(neigh[v], deg[v], u) != NULL);
        } else {
            return (fast_search(neigh[u], deg[u], v) != NULL);
        }
    }
    // Backup graph [sizeof(igraph_integer_t) bytes per edge]
    igraph_integer_t* backup(igraph_integer_t *here = NULL);
    // Restore from backup. Assume that degrees haven't changed
    void restore(igraph_integer_t* back);
    // Resplace with hard backup.
    void replace(igraph_integer_t* _hardbackup);
    // Backup degs of graph
    igraph_integer_t* backup_degs(igraph_integer_t *here = NULL);
    // Restore degs from neigh[]. Need last degree, though
    void restore_degs(igraph_integer_t last_degree);
    // Restore degs[] from backup. Assume that links[] has only been permuted
    void restore_degs_only(igraph_integer_t* backup_degs);
    // Restore degs[] and neigh[]. Assume that links[] has only been permuted
    void restore_degs_and_neigh(igraph_integer_t* backup_degs);
    // sort adjacency lists
    void sort();
    // count cycles passing through vertex v
    //int cycles(int v);
    // remove vertex (i.e. remove all edges adjacent to vertex)
    //void remove_vertex(int v);

    // For debug purposes : verify validity of the graph (symetry, simplicity)
#define VERIFY_NORMAL  0
#define VERIFY_NONEIGH 1
#define VERIFY_NOARCS  2
    bool verify(int mode = VERIFY_NORMAL);

    /*___________________________________________________________________________________
      Not to use anymore : use graph_molloy_hash class instead


    public:
      // Shuffle. returns number of swaps done.
      void shuffle(long);
      // Connected Shuffle
      long connected_shuffle(long);
      // Get caracteristic K
      double eval_K(int quality = 100);
      // Get effective K
      double effective_K(int K, int quality = 10000);
      // Test window
      double window(int K, double ratio);
      // Try to shuffle n times. Return true if at the end, the graph was still connected.
      bool try_shuffle(int T, int K);

    //___________________________________________________________________________________
    */

    /*___________________________________________________________________________________
      Not to use anymore : replaced by vertex_betweenness()     22/04/2005

      // shortest paths where vertex is an extremity
      long long *vertex_betweenness_usp(bool trivial_path);
      // shortest paths where vertex is an extremity
      long long *vertex_betweenness_rsp(bool trivial_path);
      // same, but when multiple shortest path are possible, average the weights.
      double *vertex_betweenness_asp(bool trivial_path);
    //___________________________________________________________________________________
    */

};

} // namespace gengraph

#endif //GRAPH_MOLLOY_OPT_H
