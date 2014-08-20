/* 
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

/**
 * File: bisbm_communities.c
 * -------------------------
 * Author: Robbie Ostrow <ostrowr@gmail.com>
 *   (Written while working at Lab41/In-Q-Tel)
 * Date: August 2014 
 *
 * This file provides the implementation of a bipartite stochastic block model 
 * (biSBM) algorithm for community detection in bipartite networks.
 *
 * The algorithm was originally described by Daniel Larremore, Aaron Clauset,
 * and Abigail Jacobs: Efficiently inferring community structure in bipartite 
 * networks.
 *
 * Daniel B. Larremore, Aaron Clauset, and Abigail Z. Jacobs, 
 *    Physical Review E 90(1), 012805 (2014).
 *
 * A free version of the paper can be found at http://arxiv.org/abs/1403.2933
 */


#include "igraph_interface.h"
#include "igraph_community.h"
#include "igraph_random.h"
#include "igraph_bipartite.h"
#include "igraph_conversion.h"
#include "igraph_error.h"
#include "igraph_statusbar.h"
#include "config.h"

#include <limits.h>


#define true 1
#define false 0

#define INTERCOMMUNITY(hk, c_a, c_b) (MATRIX(*(hk->inter_comm_edges), c_a, c_b - hk->a))


/**
 * TODO
 *     * possibly replace c arrays with igraph_vector_t/igraph_vector_integer_t
 *     * implement non degree-corrected version
 *     * do something to deal with degree 0 vertices (put them all in community -1 or something?)
 *     * possibly modify the algorithm to a weight corrected version
 *     * parallelize on each vertex
 *     * decompose Housekeeping struct
 */


/**
 * Struct: Housekeeping
 * --------------------
 * This struct contains all information needed to calculate the scores
 * of a partition.
 */
typedef struct {
  // number of communities of type a, number of communities of type b, number of vertices
  int a, b, size; 
  // community membership list
  int *partition; 
  // the graph in adjacency list format
  igraph_vector_t **adj_list; 
  // For each vertex i, false if types[i] is in group a, true if group b
  igraph_vector_bool_t *types; 
  // a matrix in which cell ab is the number of edges between community a and b
  igraph_matrix_t *inter_comm_edges;
  // For each community, the sum of the degree of all vertices in the community.
  igraph_real_t *comm_tot_degree; 
} Housekeeping;


/**
 * Struct: Swaprecord
 * ------------------
 * A Swaprecord contains the vertex id (v), its original community (src)
 * and the destination community (dest).
 */
typedef struct {
  int v, src, dest;
} Swaprecord;


/**
 * Function: score_partition
 * -------------------------
 * Given a partition stored in hk, scores the partition using the 
 * log-likelihood metric described in equation (16) of the paper. 
 *   Parameters:
 *     Housekeeping *hk: A struct Housekeeping containing all graph 
 *                       and community information
 *   Returns: 
 *     The score of the partition stored in hk.
 * 
 * TODO: implement non-degree corrected (equation 9).
 * TODO: implement change in score instead of total score for an optimization.
 * 
 */
static igraph_real_t score_partition(Housekeeping *hk){
  igraph_real_t score = 0;
  for (int c_a = 0; c_a < hk->a; c_a++){
    for (int c_b = hk->a; c_b < hk->a + hk->b; c_b++){
      igraph_real_t inter_community = INTERCOMMUNITY(hk, c_a, c_b);
      igraph_real_t c_a_sum = hk->comm_tot_degree[c_a];
      igraph_real_t c_b_sum = hk->comm_tot_degree[c_b];
      if (!(int) inter_community || !(int) c_a_sum || !(int) c_b_sum) return -INFINITY;
      // this line takes half the total time... possibly implement a cache?
      score += inter_community * log(inter_community / (c_a_sum * c_b_sum)); 
    }
  }
  return score;
}


/**
 * Function: make_swap
 * -------------------
 * Swaps vertex v from its current community to the new community  `to`.
 * All changes are reflected in hk.
 *   Parameters:
 *     Housekeeping *hk: A struct Housekeeping containing all graph 
 *                       and community information
 *     int v: The index of the vertex to be swapped.
 *     int to: the community to which v is to be swapped.
 */
static void make_swap(Housekeeping *hk, int v, int to){
  int from = hk->partition[v];
  if (to == from) return;
  hk->partition[v] = to;
  igraph_vector_t *neighbors = hk->adj_list[v];
  int degree = igraph_vector_size(neighbors);
  hk->comm_tot_degree[from] -= degree;
  hk->comm_tot_degree[to] += degree;
  int c_1 = from;
  for (int neigh = 0; neigh < degree; neigh++){
    int neighbor_comm = hk->partition[(int) VECTOR(*neighbors)[neigh]];
    int c_2 = neighbor_comm;
    if (c_1 < c_2){
      INTERCOMMUNITY(hk, c_1, c_2) -= 1;
      INTERCOMMUNITY(hk, to, c_2) += 1;
    }
    else{
      INTERCOMMUNITY(hk, c_2, c_1) -= 1;
      INTERCOMMUNITY(hk, c_2, to) += 1;
    }
  }
}


/**
 * Function: score_swap
 * --------------------
 * A naive scoring method that actually makes the swap that we
 * aim to score, scores the partition, then swaps back. Returns
 * the score.
 *   Parameters:
 *     Housekeeping *hk: A struct Housekeeping containing all graph 
 *                       and community information
 *     int v: The index of the vertex to be swapped for scoring.
 *     int to: the community to which v is to be swapped for scoring.
 *   Returns: 
 *     The score of the partition after the swap specified by the parameters.
 *  
 * TODO: optimize
 */
static double score_swap(Housekeeping *hk, int v, int to){
  int tmp = hk->partition[v];
  make_swap(hk, v, to);
  double score = score_partition(hk);
  make_swap(hk, v, tmp);
  return score;
}


/**
 * Function: score_swaps
 * ---------------------
 * Given a vertex v, populates dest with the best community to change to.
 * Returns the score if swapped to dest.
 *   Parameters:
 *     Housekeeping *hk: A struct Housekeeping containing all graph 
 *                       and community information
 *     int v: The index of the vertex to be swapped for scoring.
 *     int *dest: A pointer to an uninitialized int where the best swap will be stored.
 *   Returns: 
 *     The score of the optimal swap for v. 
 *  
 */
static double score_swaps(Housekeeping *hk, int v, int *dest){
  int from = hk->partition[v];
  int begin = 0;
  int end = hk->a;
  if (from >= hk->a){
    begin = hk->a; 
    end = hk->a + hk->b;
  }
  double best_score = -INFINITY;
  for (int to = begin; to < end; to++){
    if (to == from) continue;
    double new_score = score_swap(hk, v, to);
    if (new_score > best_score){
      best_score = new_score;
      *dest = to;
    }
  }
  return best_score;
}


/**
 * Function: make_best_swap
 * ------------------------
 * This function tries all possible swaps of all unused vertices to all communities,
 * then makes the best one. Returns the score after this optimal swap.
 * Takes a housekeeping struct, a bool array representing whether each vertex has been 
 * used, an array of swaprecords to give the ability to rewind to the best swap, and 
 * and int i representing the swap number that we're on. 
 *   Parameters:
 *     Housekeeping *hk: A struct Housekeeping containing all graph 
 *                       and community information
 *     bool *used: An array of length hk->size where used[i] is false if vertex i is in 
 *                 group a, true if group b.
 *     Swaprecord *swaprecord: An array of swaprecords in order, each representing one swap.
 *     int i: the number of vertices already swapped, and the swaprecord index.
 *   Returns:
 *     The score of the partition after the best swap.
 *  
 */
static double make_best_swap(Housekeeping *hk, igraph_bool_t *used, Swaprecord *swaprecord, int i){
  int max_v = -1;
  double max_score = -INFINITY;
  int max_swap = -1;
  for (int v = 0; v < hk->size; v++){
    if (used[v]) continue;
    int dest;
    double score = score_swaps(hk, v, &dest);
    if (score > max_score){      
      max_score = score;
      max_swap = dest;
      max_v = v;
    }
  }
  if (max_v == -1){
    // don't make any swaps -- we assume an empty community is illegal.
    // a little hackish here -- we just record an empty swap on vertex 0.
    swaprecord[i].v = 0;
    swaprecord[i].src = hk->partition[0];
    swaprecord[i].dest = hk->partition[0];
    return max_score;
  }
  swaprecord[i].v = max_v;
  swaprecord[i].src = hk->partition[max_v];
  swaprecord[i].dest = max_swap;
  make_swap(hk, max_v, max_swap);
  used[max_v] = true;
  return max_score;
}


/**
 * Function: rewind_swaps
 * ----------------------
 * Given an array of Swaprecords along with the index in that array best_swap,
 * swaps vertices back and updates housekeeping to reflect the best swap for the
 * current iteration.
 */
static void rewind_swaps(Housekeeping *hk, Swaprecord *swaprecord, int best_swap){
  for (int i = hk->size - 1; i > best_swap; i--){
    make_swap(hk, swaprecord[i].v, swaprecord[i].src);
  }
}


/**
 * Function: run_iteration
 * -----------------------
 * Runs the algorithm one iteration forward. If it improves the score, 
 * return false, otherwise return true.
 */
static igraph_bool_t run_iteration(Housekeeping *hk, double init_score){
  igraph_bool_t used[hk->size];
  Swaprecord swaprecord[hk->size];
  int best_swap = -1;
  for (int i = 0; i < hk->size; i++) used[i] = false;
  for (int i = 0; i < hk->size; i++){ 
    double new_score = make_best_swap(hk, used, swaprecord, i);
    if (new_score > init_score){
      init_score = new_score;
      best_swap = i;
    }
  }
  rewind_swaps(hk, swaprecord, best_swap);
  if (best_swap == -1)
    return true; // TODO make clearer
  return false;
}


/**
 * Function: run_algorithm
 * -----------------------
 * Run the bipartite SBM algorithm a maximum of max_iters iterations.
 */
static double run_algorithm(Housekeeping *hk, int max_iters){
  double score = score_partition(hk);
  if (!max_iters) max_iters = INT_MAX;
  for (int i = 0; i < max_iters; i++){
    // TODO: check that this statusf is doing something.
    igraph_statusf("Beginning iteration %d. Current score: %f.\n", 0, i + 1, score);
    igraph_bool_t is_last_iteration = run_iteration(hk, score);
    score = score_partition(hk);
    if (is_last_iteration){
      igraph_statusf("Score has not improved. Final score: %f\n", 0, score);
      return score;
    }
  }
  return score;
}


/**
 * Function: initialize_types
 * --------------------------
 * Attempts to find a bipartite mapping of the graph, and loads such
 * a mapping (if it exists) into hk->types. 
 */
static int initialize_types(Housekeeping *hk, igraph_t *graph){
  IGRAPH_VECTOR_BOOL_INIT_FINALLY(hk->types, hk->size);
  igraph_bool_t is_bipartite;
  IGRAPH_CHECK(igraph_is_bipartite(graph, &is_bipartite, hk->types));
  if (!is_bipartite){
    IGRAPH_ERROR("Input graph is not bipartite", IGRAPH_EINVAL);
  }
  return 0;
}


/**
 * Function: initialize_partition
 * ------------------------------
 * Randomly initializes hk->partition such that all vertices of type
 * a are assigned a community 0 < c < k_a and all vertices of type b
 * are assigned a community k_a <= c < k_b.
 */
static int initialize_partition(Housekeeping *hk){
  int *partition = malloc(sizeof(int) * hk->size);
  if (partition==NULL) IGRAPH_ERROR("Malloc failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(free, partition);
  igraph_rng_t *rng = igraph_rng_default();
  // assign seed to some constant for repeatable results.
  int seed = time(NULL); 
  IGRAPH_CHECK(igraph_rng_seed(rng, seed));
  for (int i = 0; i < hk->size; i++){
    if (!VECTOR(*(hk->types))[i]) // type 0
      partition[i] = igraph_rng_get_integer(rng, 0, hk->a - 1);
    else // type 1
      partition[i] = igraph_rng_get_integer(rng, hk->a, hk->a + hk->b - 1);
  }
  hk->partition = partition;
  return 0;
}


/**
 * Function: initialize_neighbors
 * ---------------------------
 * Assigns hk->adj_list as an adjacency list representation of the graph.
 *
 * TODO: find out if I can access neighbors as quickly without this data
 * structure.
 */
static int initialize_neighbors(Housekeeping *hk, igraph_t *graph){
  igraph_vector_t **neighbors = malloc(sizeof(igraph_vector_t *) * hk->size);
  if (neighbors==NULL) IGRAPH_ERROR("Malloc failed", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(free, neighbors);
  for (int v = 0; v < hk->size; v++){
    neighbors[v] = malloc(sizeof(igraph_vector_t));
    if (neighbors[v]==NULL) IGRAPH_ERROR("Malloc failed", IGRAPH_ENOMEM);
    IGRAPH_VECTOR_INIT_FINALLY(neighbors[v], 0);
    //IGRAPH_CHECK(igraph_vector_init(neighbors[v], 0));
    IGRAPH_CHECK(igraph_neighbors(graph, neighbors[v], v, IGRAPH_ALL));
  }
  hk->adj_list = neighbors;
  return 0;
}


/**
 * Function: initialize_inter_comm
 * -------------------------------
 * Initializes an A x B matrix hk->inter_comm_edges such that 
 * inter_comm_edges[A][B] = the number of intercommunity edges from A to B.
 * since the graph is undirected, only initializes the upper trianglular portion. 
 */
static int initialize_inter_comm(Housekeeping *hk, igraph_t *graph){
  igraph_matrix_t mat;
  IGRAPH_MATRIX_INIT_FINALLY(&mat, hk->size, hk->size);
  IGRAPH_CHECK(igraph_get_adjacency(graph, &mat, IGRAPH_GET_ADJACENCY_UPPER, false));
  IGRAPH_MATRIX_INIT_FINALLY(hk->inter_comm_edges, hk->a, hk->b);
  igraph_matrix_null(hk->inter_comm_edges);
  for (int row = 0; row < hk->size; row++){
    for (int col = row + 1; col < hk->size; col++){
      if (MATRIX(mat, row, col)){
        int group_r = hk->partition[row];
        int group_c = hk->partition[col];
        if (group_r < hk->a && group_c >= hk->a)
          INTERCOMMUNITY(hk, group_r, group_c)++;
        else if (group_r >= hk->a && group_r < hk->a)
          INTERCOMMUNITY(hk, group_c, group_r)++;
      }
    }
  }
  igraph_matrix_destroy(&mat);
  IGRAPH_FINALLY_CLEAN(1); // TODO: make sure the stack is behaving properly
  return 0;
}


/**
 * Function: initialize_degree_sums
 * --------------------------------
 * Initializes an array hk->comm_tot_degree where comm_tot_degree[i] = the total degree of all vertices
 * in community i.
 */
static int initialize_degree_sums(Housekeeping *hk){
  igraph_real_t *comm_degree = calloc(hk->a + hk->b, sizeof(igraph_real_t));
  if (comm_degree==NULL) IGRAPH_ERROR("Calloc failed", IGRAPH_ENOMEM);
  for (int v = 0; v < hk->size; v++){
    comm_degree[hk->partition[v]] += igraph_vector_size(hk->adj_list[v]);
  }
  hk->comm_tot_degree = comm_degree;
  return 0;
}


/**
 * Function: initialize_housekeeping
 * ---------------------------------
 * Initializes all data structures needed to run the algorithm.
 * TODO: allow the user to supply an id or something in community a, so they aren't assigned arbitrarily.
 */
static int initialize_housekeeping(Housekeeping *hk, igraph_t *graph, igraph_integer_t k_a, igraph_integer_t k_b){
  hk->a = k_a;
  hk->b = k_b;
  hk->size = igraph_vcount(graph);
  IGRAPH_CHECK(initialize_types(hk, graph));
  IGRAPH_CHECK(initialize_partition(hk));
  IGRAPH_CHECK(initialize_neighbors(hk, graph));
  IGRAPH_CHECK(initialize_inter_comm(hk, graph));
  IGRAPH_CHECK(initialize_degree_sums(hk));
  return 0;
}


/**
 * Function: free_housekeeping
 * ---------------------------
 * Frees all allocated memory in housekeeping.
 */
static void free_housekeeping(Housekeeping *hk){
  int freed = 0;
  free(hk->partition);
  freed++;
  for (int i = 0; i < hk->size; i++){
    igraph_vector_destroy(hk->adj_list[i]);
    free(hk->adj_list[i]);
    freed += 2;
  }
  free(hk->adj_list);
  freed++;
  igraph_vector_bool_destroy(hk->types);
  freed++;
  igraph_matrix_destroy(hk->inter_comm_edges);
  freed++;
  free(hk->comm_tot_degree);
  freed++;
  IGRAPH_FINALLY_CLEAN(freed); // could hardcode freed instead
}


/**
 * \function igraph_community_bipartite_sbm
 *
 * 
 * This function is an implementation of the degree-corrected bipartite 
 * stochastic block model community detection algorithm (biSBM) by 
 * Larremore et. al in
 * Larremore, D., Clauset, A. and Jacobs, A. (2014). Efficiently 
 *    inferring community structure in bipartite networks. 
 *    Phys. Rev. E, [online] 90, p.012805. 
 *
 * A free version of the paper can be read here: http://arxiv.org/abs/1403.2933
 *
 * </para><para>
 * Currently, only the degree-corrected version of their algorithm is
 * implemented. 
 *
 * \param graph The input graph, assumed to be bipartite and undirected.
 * \param membership Pointer to an initialized igraph_vector_t, the 
 *    community memberships will be stored here.
 * \param k_a The number of communities of the first type.
 * \param k_b The number of communities of the second type. Note that you
 *    cannot currently specify which group is which.
 * \param max_iters The maximum number of iterations to run the algorithm. If 
 * 0, then the algorithm will run until it converges.
 * \return The maximum score attained. Since this algorithm is stochastic,
 * it can return different scores on different runs. The maximum score among
 * all of those runs is the optimal community structure.
 *
 * Time complexity: roughly O(N_aK_a(K_a + <k>)) + O(N_bK_b(K_b + <k>)) where 
 * N_a is the number of vertices of type a, K_a is the number of communities of 
 * type a, and <k> is the average degree of each vertex. 
 *
 * \example examples/simple/igraph_community_bipartite_sbm.c
 */
int igraph_community_bipartite_sbm(igraph_t *graph, igraph_vector_t *membership, 
                                   igraph_integer_t k_a, igraph_integer_t k_b, 
                                   igraph_integer_t max_iters){
  Housekeeping hk;
  igraph_vector_bool_t types;
  hk.types = &types;
  igraph_matrix_t inter_comm_edges;
  hk.inter_comm_edges = &inter_comm_edges;
  IGRAPH_CHECK(initialize_housekeeping(&hk, graph, k_a, k_b));
  run_algorithm(&hk, max_iters);
  if (igraph_vector_size(membership) != hk.size)
    IGRAPH_CHECK(igraph_vector_resize(membership, hk.size));
  for (int i = 0; i < hk.size; i++)
    VECTOR(*membership)[i] = hk.partition[i];
  free_housekeeping(&hk);
  return 0;
}