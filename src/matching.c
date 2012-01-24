/* -*- mode: C -*-  */
/* vim:set ts=2 sts=2 sw=2 et: */
/* 
   IGraph library.
   Copyright (C) 2007-2011  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#include "config.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_matching.h"
#include "igraph_structural.h"

/* #define MATCHING_DEBUG */

#ifdef _MSC_VER
/* MSVC does not support variadic macros */
#include <stdarg.h>
static void debug(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
#ifdef MATCHING_DEBUG
  vfprintf(stderr, fmt, args);
#endif
  va_end(args);
}
#else
#  ifdef MATCHING_DEBUG
#    define debug(...) fprintf(stderr, __VA_ARGS__)
#  else
#    define debug(...) 
#  endif
#endif

/**
 * \function igraph_is_matching
 * Checks whether the given matching is valid for the given graph.
 *
 * This function checks a matching vector and verifies whether its length
 * matches the number of vertices in the given graph, its values are between
 * -1 (inclusive) and the number of vertices (exclusive), and whether there
 * exists a corresponding edge in the graph for every matched vertex pair.
 * For bipartite graphs, it also verifies whether the matched vertices are
 * in different parts of the graph.
 *
 * \param graph The input graph. It can be directed but the edge directions
 *              will be ignored.
 * \param types If the graph is bipartite and you are interested in bipartite
 *              matchings only, pass the vertex types here. If the graph is
 *              non-bipartite, simply pass \c NULL.
 * \param matching The matching itself. It must be a vector where element i
 *                 contains the ID of the vertex that vertex i is matched to,
 *                 or -1 if vertex i is unmatched.
 * \param result Pointer to a boolean variable, the result will be returned
 *               here.
 *
 * \sa \ref igraph_is_maximal_matching() if you are also interested in whether
 *     the matching is maximal (i.e. non-extendable).
 *
 * Time complexity: O(|V|+|E|) where |V| is the number of vertices and
 * |E| is the number of edges.
 * 
 * \example examples/simple/igraph_maximum_bipartite_matching.c
 */
int igraph_is_matching(const igraph_t* graph,
    const igraph_vector_bool_t* types, const igraph_vector_long_t* matching,
    igraph_bool_t* result) {
  long int i, j, no_of_nodes = igraph_vcount(graph);
  igraph_bool_t conn;

  /* Checking match vector length */
  if (igraph_vector_long_size(matching) != no_of_nodes) {
    *result = 0; return IGRAPH_SUCCESS;
  }

  for (i = 0; i < no_of_nodes; i++) {
    j = VECTOR(*matching)[i];

    /* Checking range of each element in the match vector */
    if (j < -1 || j >= no_of_nodes) {
      *result = 0; return IGRAPH_SUCCESS;
    }
    /* When i is unmatched, we're done */
    if (j == -1)
      continue;
    /* Matches must be mutual */
    if (VECTOR(*matching)[j] != i) {
      *result = 0; return IGRAPH_SUCCESS;
    }
    /* Matched vertices must be connected */
    IGRAPH_CHECK(igraph_are_connected(graph, i, j, &conn));
    if (!conn) {
      /* Try the other direction -- for directed graphs */
      IGRAPH_CHECK(igraph_are_connected(graph, j, i, &conn));
      if (!conn) {
        *result = 0; return IGRAPH_SUCCESS;
      }
    }
  }

  if (types != 0) {
    /* Matched vertices must be of different types */
    for (i = 0; i < no_of_nodes; i++) {
      j = VECTOR(*matching)[i];
      if (j == -1)
        continue;
      if (VECTOR(*types)[i] == VECTOR(*types)[j]) {
        *result = 0; return IGRAPH_SUCCESS;
      }
    }
  }

  *result = 1;
  return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_maximal_matching
 * Checks whether a matching in a graph is maximal.
 *
 * A matching is maximal if and only if there exists no unmatched vertex in a
 * graph such that one of its neighbors is also unmatched.
 *
 * \param graph The input graph. It can be directed but the edge directions
 *              will be ignored.
 * \param types If the graph is bipartite and you are interested in bipartite
 *              matchings only, pass the vertex types here. If the graph is
 *              non-bipartite, simply pass \c NULL.
 * \param matching The matching itself. It must be a vector where element i
 *                 contains the ID of the vertex that vertex i is matched to,
 *                 or -1 if vertex i is unmatched.
 * \param result Pointer to a boolean variable, the result will be returned
 *               here.
 *
 * \sa \ref igraph_is_matching() if you are only interested in whether a
 *     matching vector is valid for a given graph.
 *
 * Time complexity: O(|V|+|E|) where |V| is the number of vertices and
 * |E| is the number of edges.
 * 
 * \example examples/simple/igraph_maximum_bipartite_matching.c
 */
int igraph_is_maximal_matching(const igraph_t* graph,
    const igraph_vector_bool_t* types, const igraph_vector_long_t* matching,
    igraph_bool_t* result) {
  long int i, j, n, no_of_nodes = igraph_vcount(graph);
  igraph_vector_t neis;
  igraph_bool_t valid;

  IGRAPH_CHECK(igraph_is_matching(graph, types, matching, &valid));
  if (!valid) {
    *result = 0; return IGRAPH_SUCCESS;
  }

  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  valid = 1;
  for (i = 0; i < no_of_nodes; i++) {
    j = VECTOR(*matching)[i];
    if (j != -1)
      continue;

    IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_ALL));
    n = igraph_vector_size(&neis);
    for (j = 0; j < n; j++) {
      if (VECTOR(*matching)[(long int)VECTOR(neis)[j]] == -1) {
        if (types == 0 ||
            VECTOR(*types)[i] != VECTOR(*types)[(long int)VECTOR(neis)[j]]) {
          valid = 0; break;
        }
      }
    }
  }

  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(1);

  *result = valid;
  return IGRAPH_SUCCESS;
}

int igraph_i_maximum_bipartite_matching_unweighted(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_vector_long_t* matching);
int igraph_i_maximum_bipartite_matching_weighted(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights);

#define MATCHED(v) (VECTOR(match)[v] != -1)
#define UNMATCHED(v) (!MATCHED(v))

/**
 * \function igraph_maximum_bipartite_matching
 * Calculates a maximum matching in a bipartite graph.
 *
 * A matching in a bipartite graph is a partial assignment of vertices
 * of the first kind to vertices of the second kind such that each vertex of
 * the first kind is matched to at most one vertex of the second kind and
 * vice versa, and matched vertices must be connected by an edge in the graph.
 * The size (or cardinality) of a matching is the number of edges.
 * A matching is a maximum matching if there exists no other matching with
 * larger cardinality. For weighted graphs, a maximum matching is a matching
 * whose edges have the largest possible total weight among all possible
 * matchings.
 *
 * </para><para>
 * Maximum matchings in bipartite graphs are found by the push-relabel algorithm
 * with greedy initialization and a global relabeling after every n/2 steps where
 * n is the number of vertices in the graph.
 *
 * </para><para>
 * References: Cherkassky BV, Goldberg AV, Martin P, Setubal JC and Stolfi J:
 * Augment or push: A computational study of bipartite matching and
 * unit-capacity flow algorithms. ACM Journal of Experimental Algorithmics 3,
 * 1998.
 *
 * </para><para>
 * Kaya K, Langguth J, Manne F and Ucar B: Experiments on push-relabel-based
 * maximum cardinality matching algorithms for bipartite graphs. Technical
 * Report TR/PA/11/33 of the Centre Europeen de Recherche et de Formation
 * Avancee en Calcul Scientifique, 2011.
 *
 * \param graph The input graph. It can be directed but the edge directions
 *              will be ignored.
 * \param types Boolean vector giving the vertex types of the graph.
 * \param matching_size The size of the matching (i.e. the number of matched
 *                      vertex pairs will be returned here). It may be \c NULL
 *                      if you don't need this.
 * \param matching_weight The weight of the matching if the edges are weighted,
 *                        or the size of the matching again if the edges are
 *                        unweighted. It may be \c NULL if you don't need this.
 * \param matching The matching itself. It must be a vector where element i
 *                 contains the ID of the vertex that vertex i is matched to,
 *                 or -1 if vertex i is unmatched.
 * \param weights A null pointer (=no edge weights), or a vector giving the
 *                weights of the edges.
 * \return Error code.
 *
 * Time complexity: O(sqrt(|V|) |E|) according to the technical report
 * referenced above.
 * 
 * \example examples/simple/igraph_maximum_bipartite_matching.c
 */
int igraph_maximum_bipartite_matching(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights) {
  if (weights == 0) {
    IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted(graph, types,
        matching_size, matching));
    if (matching_weight != 0) {
      *matching_weight = *matching_size;
    }
    return IGRAPH_SUCCESS;
  } else {
    return igraph_i_maximum_bipartite_matching_weighted(graph, types,
        matching_size, matching_weight, matching, weights);
  }
}

int igraph_i_maximum_bipartite_matching_unweighted_relabel(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_vector_t* labels,
    igraph_vector_long_t* matching, igraph_bool_t smaller_set);

/**
 * Finding maximum bipartite matchings on bipartite graphs using the
 * push-relabel algorithm.
 *
 * The implementation follows the pseudocode in Algorithm 1 of the
 * following paper:
 * 
 * Kaya K, Langguth J, Manne F and Ucar B: Experiments on push-relabel-based
 * maximum cardinality matching algorithms for bipartite graphs. Technical
 * Report TR/PA/11/33 of CERFACS (Centre Européen de Recherche et de Formation
 * Avancée en Calcul Scientifique).
 * http://www.cerfacs.fr/algor/reports/2011/TR_PA_11_33.pdf
 */
int igraph_i_maximum_bipartite_matching_unweighted(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_vector_long_t* matching) {
  long int i, j, k, n, no_of_nodes = igraph_vcount(graph);
  long int num_matched;             /* number of matched vertex pairs */
  igraph_vector_long_t match;       /* will store the matching */
  igraph_vector_t labels;           /* will store the labels */
  igraph_vector_t neis;             /* used to retrieve the neighbors of a node */
  igraph_dqueue_long_t q;           /* a FIFO for push ordering */
  igraph_bool_t smaller_set;        /* denotes which part of the bipartite graph is smaller */
  long int size_of_smaller_set;
  long int label_changed = 0;       /* Counter to decide when to run a global relabeling */
  long int relabeling_freq = no_of_nodes / 2;

  /* We will use:
   * - FIFO push ordering
   * - global relabeling frequency: n/2 steps where n is the number of nodes
   * - simple greedy matching for initialization
   */

  /* (1) Initialize data structures */
  IGRAPH_CHECK(igraph_vector_long_init(&match, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &match);
  IGRAPH_VECTOR_INIT_FINALLY(&labels, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
  IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);

  /* (2) Initially, every node is unmatched */
  igraph_vector_long_fill(&match, -1);

  /* (3) Find an initial matching in a greedy manner.
   *     At the same time, find which side of the graph is smaller. */
  num_matched = 0; j = 0;
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i])
      j++;
    if (MATCHED(i))
      continue;
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_ALL));
    n = igraph_vector_size(&neis);
    for (j = 0; j < n; j++) {
      k = VECTOR(neis)[j];
      if (UNMATCHED(k)) {
        /* We match vertex i to vertex VECTOR(neis)[j] */
        VECTOR(match)[k] = i;
        VECTOR(match)[i] = k;
        num_matched++;
        break;
      }
    }
  }
  smaller_set = (j <= no_of_nodes/2);
  size_of_smaller_set = (smaller_set ? j : no_of_nodes-j);

  /* (4) Set the initial labeling -- lines 1 and 2 in the tech report */
  IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted_relabel(
      graph, types, &labels, &match, smaller_set));

  /*
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i] == smaller_set) {
      VECTOR(labels)[i] = 1;
    } else {
      VECTOR(labels)[i] = 0;
    }
  }
  */

  /* Smaller set now has label 1, larger set has label 0. The tech report
   * calls the smaller set "columns" and the larger set "rows". */

  /* (5) Fill the push queue with the unmatched nodes from the smaller set. */
  for (i = 0; i < no_of_nodes; i++) {
    if (UNMATCHED(i) && VECTOR(*types)[i] == smaller_set)
      IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
  }

  /* (6) Main loop from the referenced tech report -- lines 4--13 */
  label_changed = 0;
  while (!igraph_dqueue_long_empty(&q)) {
    long int v = igraph_dqueue_long_pop(&q);             /* Line 13 */
    long int u = -1, label_u = 2 * no_of_nodes;
    long int w;

    if (label_changed >= relabeling_freq) {
      /* Run global relabeling */
      IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted_relabel(
            graph, types, &labels, &match, smaller_set));
      label_changed = 0;
    }

    debug("Considering vertex %ld\n", v);

    /* Line 5: find row u among the neighbors of v s.t. label(u) is minimal */
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, IGRAPH_ALL));
    n = igraph_vector_size(&neis);
    for (i = 0; i < n; i++) {
      if (VECTOR(labels)[(long int)VECTOR(neis)[i]] < label_u) {
        u = VECTOR(neis)[i];
        label_u = VECTOR(labels)[u];
        label_changed++;
      }
    }

    debug("  Neighbor with smallest label: %ld (label=%ld)\n", u, label_u);

    if (label_u < no_of_nodes) {                         /* Line 6 */
      VECTOR(labels)[v] = VECTOR(labels)[u] + 1;         /* Line 7 */
      if (MATCHED(u)) {                                  /* Line 8 */
        w = VECTOR(match)[u];
        debug("  Vertex %ld is matched to %ld, performing a double push\n", u, w);
        if (w != v) {
          VECTOR(match)[u] = -1; VECTOR(match)[w] = -1;  /* Line 9 */
          IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));  /* Line 10 */
          debug("  Unmatching & activating vertex %ld\n", w);
          num_matched--;
        }
      }
      VECTOR(match)[u] = v; VECTOR(match)[v] = u;      /* Line 11 */
      num_matched++;
      VECTOR(labels)[u] += 2;                          /* Line 12 */
      label_changed++;
    }
  }

  /* Fill the output parameters */
  if (matching != 0) {
    IGRAPH_CHECK(igraph_vector_long_update(matching, &match));
  }
  if (matching_size != 0) {
    *matching_size = num_matched;
  }

  /* Release everything */
  igraph_dqueue_long_destroy(&q);
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&labels);
  igraph_vector_long_destroy(&match);
  IGRAPH_FINALLY_CLEAN(4);

  return IGRAPH_SUCCESS;
}

int igraph_i_maximum_bipartite_matching_unweighted_relabel(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_vector_t* labels,
    igraph_vector_long_t* match, igraph_bool_t smaller_set) {
  long int i, j, n, no_of_nodes = igraph_vcount(graph), matched_to;
  igraph_dqueue_long_t q;
  igraph_vector_t neis;

  debug("Running global relabeling.\n");

  /* Set all the labels to no_of_nodes first */
  igraph_vector_fill(labels, no_of_nodes);

  /* Allocate vector for neighbors */
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  /* Create a FIFO for the BFS and initialize it with the unmatched rows
   * (i.e. members of the larger set) */
  IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
  IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i] != smaller_set && VECTOR(*match)[i] == -1) {
      IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
      VECTOR(*labels)[i] = 0;
    }
  }

  /* Run the BFS */
  while (!igraph_dqueue_long_empty(&q)) {
    long int v = igraph_dqueue_long_pop(&q);
    long int w;

    IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, IGRAPH_ALL));

    n = igraph_vector_size(&neis);
    for (j = 0; j < n; j++) {
      w = VECTOR(neis)[j];
      if (VECTOR(*labels)[w] == no_of_nodes) {
        VECTOR(*labels)[w] = VECTOR(*labels)[v] + 1;
        matched_to = VECTOR(*match)[w];
        if (matched_to != -1 && VECTOR(*labels)[matched_to] == no_of_nodes) {
          IGRAPH_CHECK(igraph_dqueue_long_push(&q, matched_to));
          VECTOR(*labels)[matched_to] = VECTOR(*labels)[w] + 1;
        }
      }
    }
  }

  igraph_dqueue_long_destroy(&q);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(2);

  return IGRAPH_SUCCESS;
}

int igraph_i_maximum_bipartite_matching_weighted(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights) {
  IGRAPH_ERROR("maximum weighted matching on bipartite graphs not implemented yet",
      IGRAPH_UNIMPLEMENTED);
}

int igraph_maximum_matching(const igraph_t* graph, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights) {
  IGRAPH_ERROR("maximum matching on general graphs not implemented yet",
      IGRAPH_UNIMPLEMENTED);
}

#ifdef MATCHING_DEBUG
#undef MATCHING_DEBUG
#endif


