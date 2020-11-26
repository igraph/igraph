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
#include "gengraph_qsort.h"

#include <assert.h>
#include "gengraph_random.h"

namespace gengraph {

// This class handles graphs with a constant degree sequence.

class graph_molloy_opt {

private:
    // Random generator
    KW_RNG::RNG rng;
    // Number of vertices
    int n;
    //Number of arcs ( = #edges * 2 )
    int a;
    // The degree sequence of the graph
    int *deg;
    // The array containing all links
    int *links;
    // The array containing pointers to adjacency list of every vertices
    int **neigh;
    // Allocate memory according to degree_sequence (for constructor use only!!)
    void alloc(degree_sequence &);
    // Compute #edges
    inline void refresh_nbarcs() {
        a = 0;
        for (int* d = deg + n; d != deg; ) {
            a += *(--d);
        }
    }
    // Build neigh with deg and links
    void compute_neigh();
    // Swap edges. The swap MUST be valid !!!
    inline void swap_edges(int from1, int to1, int from2, int to2) {
        fast_rpl(neigh[from1], to1, to2);
        fast_rpl(neigh[from2], to2, to1);
        fast_rpl(neigh[to1], from1, from2);
        fast_rpl(neigh[to2], from2, from1);
    }

    // Swap edges only if they are simple. return false if unsuccessful.
    bool swap_edges_simple(int, int, int, int);
    // Test if vertex is in an isolated component of size<K
    bool isolated(int v, int K, int *Kbuff, bool *visited);
    // Pick random edge, and gives a corresponding vertex
    inline int pick_random_vertex() {
        return links[my_random() % a];
    };
    // Pick random neighbour
    inline int* random_neighbour(const int v) {
        return neigh[v] + (my_random() % deg[v]);
    };
    // Returns complexity of isolation test
    long effective_isolated(int v, int K, int *Kbuff, bool *visited);
    // Depth-Exploration. Returns number of steps done. Stops when encounter vertex of degree > dmax.
    void depth_isolated(int v, long &calls, int &left_to_explore, int dmax, int * &Kbuff, bool *visited);
    // breadth-first search. Store the distance (modulo 3)  in dist[]. Returns eplorated component size.
    int width_search(unsigned char *dist, int *buff, int v0 = 0, int toclear = -1);
    // depth-first search.
    int depth_search(bool *visited, int *buff, int v0 = 0);
    // breadth-first search that count the number of shortest paths going from src to each vertex
    int breadth_path_search(int src, int *buff, double *paths, unsigned char *dist);
    // Used by traceroute_sample() ONLY
    void add_traceroute_edge(int, int, int*, double** red = NULL, double t = 1.0);
    // Used by traceroute() and betweenness(). if newdeg[]=NULL, do not discover edges.
    // breadth_path_search() must have been called to give the corresponding buff[],dist[],paths[] and nb_vertices
    void explore_usp(double *target, int nb_vertices, int *buff, double *paths, unsigned char *dist, int *newdeg = NULL, double **edge_redudancy = NULL);
    void explore_asp(double *target, int nb_vertices, int *buff, double *paths, unsigned char *dist, int *newdeg = NULL, double **edge_redudancy = NULL);
    void explore_rsp(double *target, int nb_vertices, int *buff, double *paths, unsigned char *dist, int *newdeg = NULL, double **edge_redudancy = NULL);
    // Return component indexes where vertices belong to, starting from 0,
    // sorted by size (biggest component has index 0)
    int *components(int *comp = NULL);
    // pick k random vertices of degree > 0.
    int *pick_random_vertices(int &k, int *output = NULL, int nb_v = -1, int *among = NULL);

public:
    // neigh[]
    inline int** neighbors() {
        return neigh;
    };
    // deg[]
    inline int* degrees() {
        return deg;
    };
    //adjacency list of v
    inline int* operator[](const int v) {
        return neigh[v];
    };
    //degree of v
    inline int degree(const int v) {
        return deg[v];
    };
    //compare adjacency lists
    inline int compare(const int v, const int w) {
        return deg[v] == deg[w] ? lex_comp(neigh[v], neigh[w], deg[v]) : (deg[v] > deg[w] ? -1 : 1);
    };
    // Detach deg[] and neigh[]
    void detach();
    // Destroy deg and links
    ~graph_molloy_opt();
    // Create graph from file (stdin not supported unless rewind() possible)
    graph_molloy_opt(FILE *f);
    // Allocate memory for the graph. Create deg and links. No edge is created.
    graph_molloy_opt(degree_sequence &);
    // Create graph from hard copy
    graph_molloy_opt(int *);
    // Create hard copy of graph
    int *hard_copy();
    // Remove unused edges, updates neigh[], recreate links[]
    void clean();
    // nb arcs
    inline int nbarcs() {
        return a;
    };
    // last degree
    inline int last_degree() {
        return deg[n - 1];
    };
    // nb vertices
    inline int nbvertices() {
        return n;
    };
    // nb vertices having degree > 0
    inline int nbvertices_real() {
        int s = 0;
        for (int *d = deg + n; d-- != deg; ) if (*d) {
                s++;
            }
        return s;
    };
    // return list of vertices with degree > 0. Compute #vertices, if not given.
    int *vertices_real(int &nb_v);
    // Keep only giant component
    void giant_comp();
    // nb vertices in giant component
    int nbvertices_comp();
    // nb arcs in giant component
    int nbarcs_comp();
    // print graph in SUCC_LIST mode, in stdout
    void print(FILE *f = stdout, bool NOZERO = true);
    // Bind the graph avoiding multiple edges or self-edges (return false if fail)
    bool havelhakimi();
    // Get the graph connected  (return false if fail)
    bool make_connected();
    // Test if graph is connected
    bool is_connected();
    // Maximum degree
    int max_degree();
    // breadth-first search. Store the distance (modulo 3)  in dist[].
    void breadth_search(int *dist, int v0 = 0, int* buff = NULL);
    // is edge ?
    inline bool is_edge(const int u, const int v) {
        if (deg[v] < deg[u]) {
            return (fast_search(neigh[v], deg[v], u) != NULL);
        } else {
            return (fast_search(neigh[u], deg[u], v) != NULL);
        }
    }
    // Backup graph [sizeof(int) bytes per edge]
    int* backup(int *here = NULL);
    // Restore from backup. Assume that degrees haven't changed
    void restore(int* back);
    // Resplace with hard backup.
    void replace(int* _hardbackup);
    // Backup degs of graph
    int* backup_degs(int *here = NULL);
    // Restore degs from neigh[]. Need last degree, though
    void restore_degs(int last_degree);
    // Restore degs[] from backup. Assume that links[] has only been permuted
    void restore_degs_only(int* backup_degs);
    // Restore degs[] and neigh[]. Assume that links[] has only been permuted
    void restore_degs_and_neigh(int* backup_degs);
// WARNING : the following shuffle() algorithms are slow.
// Use graph_molloy_hash::connected_shuffle() instead.
    // "Fab" Shuffle (Optimized heuristic of Gkantsidis algo.)
    long fab_connected_shuffle(long);
    // "Optimized-Fab" Shuffle (Optimized heuristic of Gkantsidis algo, with isolated pairs)
    long opt_fab_connected_shuffle(long);
    // Gkantsidis Shuffle
    long gkantsidis_connected_shuffle(long);
    // Connected Shuffle
    long slow_connected_shuffle(long);
    // shortest paths where vertex is an extremity
    double *vertex_betweenness(int mode, bool trivial_path = false);
    // Sample the graph with traceroute-like exploration from src[] to dst[].
    // if dst[]=NULL, pick nb_dst new random destinations for each src
    double traceroute_sample(int mode, int nb_src, int *src, int nb_dst, int* dst, double *redudancy = NULL, double **edge_redudancy = NULL);
    // does one breadth-first search and returns the average_distance.
    double avg_dist(unsigned char *dist, int *buff, int v0, int &nb_vertices, int toclear = -1);
    // Number of edges needed to disconnect graph (one random instance)
    int disconnecting_edges();
    // Compute vertex covering of the graph. Warning : this modifies degs[]
    void vertex_covering();
    // Path sampling. Input is nb_dst[] and dst[]. nb_dst[v],dst[v] describe all paths (v,x)
    double path_sampling(int *nb_dst, int *dst = NULL, double *redudancies = NULL, double **edge_redudancy = NULL);
    // keep only core (tree parts are deleted). Returns number of removed vertices.
    int core();
    // try to disconnect the graph by swapping edges (with isolation tests)
    int try_disconnect(int K, int max_tries = 10000000);
    // Eric & Cun-Hui estimator
    double rho(int mode, int nb_src, int *src, int nb_dst, int *dst = NULL);
    // sort adjacency lists
    void sort();
    // sort the vertices according to their degrees (highest first) and to their adjacency lists (lexicographic)
    int* sort_vertices(int *buff = NULL);
    // count cycles passing through vertex v
    int cycles(int v);
    // remove vertex (i.e. remove all edges adjacent to vertex)
    void remove_vertex(int v);
    // pick k random vertices of degree > 0. If k \in [0,1[, k is understood as a density.
    int *pick_random_src(double k, int *nb = NULL, int* buff = NULL, int nb_v = -1, int* among = NULL);
    // pick k random vertices of degree > 0. If k \in [0,1], k is understood as a density.
    int *pick_random_dst(double k, int *nb = NULL, int* buff = NULL, int nb_v = -1, int* among = NULL);

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


