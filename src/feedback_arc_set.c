/* -*- mode: C -*-  */
/* vim:set ts=2 sts=2 sw=2 et: */
/* 
   IGraph library.
   Copyright (C) 2011  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include "igraph_centrality.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_structural.h"
#include "igraph_types.h"

int igraph_i_feedback_arc_set_undirected(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights);
int igraph_i_feedback_arc_set_eades(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights);
int igraph_i_feedback_arc_set_ip(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights);


/**
 * \ingroup structural
 * \function igraph_feedback_arc_set
 * \brief Calculates a feedback arc set of the graph using different
 *        algorithms.
 *
 * </para><para>
 * A feedback arc set is a set of edges whose removal makes the graph acyclic.
 * We are usually interested in \em minimum feedback arc sets, i.e. sets of edges
 * whose total weight is minimal among all the feedback arc sets.
 *
 * </para><para>
 * For undirected graphs, the problem is simple: one has to find a maximum weight
 * spanning tree and then remove all the edges not in the spanning tree. For directed
 * graphs, this is an NP-hard problem, and various heuristics are usually used to
 * find an approximate solution to the problem. This function implements a few of
 * these heuristics.
 *
 * \param graph  The graph object.
 * \param result An initialized vector, the result will be returned here.
 * \param weights Weight vector or NULL if no weights are specified.
 * \param algo   The algorithm to use to solve the problem if the graph is directed.
 *        Possible values:
 *        \clist
 *        \cli IGRAPH_FAS_EXACT_IP
 *          Finds a \em minimum feedback arc set using integer programming (IP).
 *          The complexity of this algorithm is exponential of course.
 *        \cli IGRAPH_FAS_APPROX_EADES
 *          Finds a feedback arc set using the heuristic of Eades, Lin and
 *          Smyth (1993). This is guaranteed to be smaller than |E|/2 - |V|/6,
 *          and it is linear in the number of edges (i.e. O(|E|)).
 *          For more details, see Eades P, Lin X and Smyth WF: A fast and effective
 *          heuristic for the feedback arc set problem. In: Proc Inf Process Lett
 *          319-323, 1993.
 *        \endclist
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL if an unknown method was specified or the weight vector
 *            is invalid.
 *
 * \example examples/simple/igraph_feedback_arc_set.c
 *
 * Time complexity: depends on \p algo, see the time complexities there.
 */
int igraph_feedback_arc_set(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights, igraph_fas_algorithm_t algo) {

  if (weights && igraph_vector_size(weights) < igraph_ecount(graph))
    IGRAPH_ERROR("cannot calculate feedback arc set, weight vector too short",
      IGRAPH_EINVAL);

  if (!igraph_is_directed(graph))
    return igraph_i_feedback_arc_set_undirected(graph, result, weights);

  switch (algo) {
    case IGRAPH_FAS_EXACT_IP:
      return igraph_i_feedback_arc_set_ip(graph, result, weights);

    case IGRAPH_FAS_APPROX_EADES:
      return igraph_i_feedback_arc_set_eades(graph, result, weights);

    default:
      IGRAPH_ERROR("Invalid algorithm", IGRAPH_EINVAL);
  }
}

/**
 * Solves the feedback arc set problem for undirected graphs.
 */
int igraph_i_feedback_arc_set_undirected(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights) {
  IGRAPH_ERROR("TODO", IGRAPH_UNIMPLEMENTED);
}

/**
 * Solves the feedback arc set problem using the heuristics of Eades et al.
 */
int igraph_i_feedback_arc_set_eades(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights) {
  long int i, j, k, v, eid, no_of_nodes=igraph_vcount(graph), nodes_left;
  igraph_dqueue_t sources, sinks;
  igraph_vector_t neis;
  igraph_vector_t indegrees, outdegrees;
  igraph_vector_t instrengths, outstrengths;
  long int* ordering;
  long int order_next_pos = 0, order_next_neg = -1;
  igraph_real_t diff, maxdiff;

  ordering = igraph_Calloc(no_of_nodes, long int);
  IGRAPH_FINALLY(igraph_free, ordering);

  IGRAPH_VECTOR_INIT_FINALLY(&indegrees, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&outdegrees, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&instrengths, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&outstrengths, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
  IGRAPH_CHECK(igraph_dqueue_init(&sinks, 0));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &sinks);

  IGRAPH_CHECK(igraph_degree(graph, &indegrees, igraph_vss_all(), IGRAPH_IN, 0));
  IGRAPH_CHECK(igraph_degree(graph, &outdegrees, igraph_vss_all(), IGRAPH_OUT, 0));

  if (weights) {
    IGRAPH_CHECK(igraph_strength(graph, &instrengths, igraph_vss_all(), IGRAPH_IN, 0, weights));
    IGRAPH_CHECK(igraph_strength(graph, &outstrengths, igraph_vss_all(), IGRAPH_OUT, 0, weights));
  } else {
    IGRAPH_CHECK(igraph_vector_update(&instrengths, &indegrees));
    IGRAPH_CHECK(igraph_vector_update(&outstrengths, &outdegrees));
  }

  /* Find initial sources and sinks */
  nodes_left = no_of_nodes;
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(indegrees)[i] == 0) {
      if (VECTOR(outdegrees)[i] == 0) {
        /* Isolated vertex, we simply ignore it */
        nodes_left--;
        ordering[i] = order_next_pos++;
      } else {
        /* This is a source */
        igraph_dqueue_push(&sources, i);
      }
    } else if (VECTOR(outdegrees)[i] == 0) {
      /* This is a sink */
      igraph_dqueue_push(&sinks, i);
    }
  }

  /* While we have any nodes left... */
  while (nodes_left > 0) {
    /* (1) Remove the sources one by one */
    while (!igraph_dqueue_empty(&sources)) {
      i=(long)igraph_dqueue_pop(&sources);
      /* Add the node to the ordering */
      ordering[i] = order_next_pos++;
      /* Exclude the node from further searches */
      VECTOR(indegrees)[i] = VECTOR(outdegrees)[i] = -1;
      /* Get the neighbors and decrease their degrees */
      IGRAPH_CHECK(igraph_incident(graph, &neis, i, IGRAPH_OUT));
      j = igraph_vector_size(&neis);
      for (i = 0; i < j; i++) {
        eid = VECTOR(neis)[i];
        k = IGRAPH_TO(graph, eid);
        if (VECTOR(indegrees)[k] <= 0) {
          /* Already removed, continue */
          continue;
        }
        VECTOR(indegrees)[k]--;
        VECTOR(instrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
        if (VECTOR(indegrees)[k] == 0)
          IGRAPH_CHECK(igraph_dqueue_push(&sources, k));
      }
      nodes_left--;
    }

    /* (2) Remove the sinks one by one */
    while (!igraph_dqueue_empty(&sinks)) {
      i=(long)igraph_dqueue_pop(&sinks);
      /* Maybe the vertex became sink and source at the same time, hence it
       * was already removed in the previous iteration. Check it. */
      if (VECTOR(indegrees)[i] < 0)
        continue;
      /* Add the node to the ordering */
      ordering[i] = order_next_neg--;
      /* Exclude the node from further searches */
      VECTOR(indegrees)[i] = VECTOR(outdegrees)[i] = -1;
      /* Get the neighbors and decrease their degrees */
      IGRAPH_CHECK(igraph_incident(graph, &neis, i, IGRAPH_IN));
      j = igraph_vector_size(&neis);
      for (i = 0; i < j; i++) {
        eid = VECTOR(neis)[i];
        k = IGRAPH_FROM(graph, eid);
        if (VECTOR(outdegrees)[k] <= 0) {
          /* Already removed, continue */
          continue;
        }
        VECTOR(outdegrees)[k]--;
        VECTOR(outstrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
        if (VECTOR(outdegrees)[k] == 0)
          IGRAPH_CHECK(igraph_dqueue_push(&sinks, k));
      }
      nodes_left--;
    }

    /* (3) No more sources or sinks. Find the node with the largest
     * difference between its out-strength and in-strength */
    v = -1; maxdiff = -IGRAPH_INFINITY;
    for (i = 0; i < no_of_nodes; i++) {
      if (VECTOR(outdegrees)[i] < 0)
        continue;
      diff = VECTOR(outstrengths)[i]-VECTOR(instrengths)[i];
      if (diff > maxdiff) {
        maxdiff = diff;
        v = i;
      }
    }
    if (v >= 0) {
      /* Remove vertex v */
      ordering[v] = order_next_pos++;
      /* Remove outgoing edges */
      IGRAPH_CHECK(igraph_incident(graph, &neis, v, IGRAPH_OUT));
      j = igraph_vector_size(&neis);
      for (i = 0; i < j; i++) {
        eid = VECTOR(neis)[i];
        k = IGRAPH_TO(graph, eid);
        if (VECTOR(indegrees)[k] <= 0) {
          /* Already removed, continue */
          continue;
        }
        VECTOR(indegrees)[k]--;
        VECTOR(instrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
        if (VECTOR(indegrees)[k] == 0)
          IGRAPH_CHECK(igraph_dqueue_push(&sources, k));
      }
      /* Remove incoming edges */
      IGRAPH_CHECK(igraph_incident(graph, &neis, v, IGRAPH_IN));
      j = igraph_vector_size(&neis);
      for (i = 0; i < j; i++) {
        eid = VECTOR(neis)[i];
        k = IGRAPH_FROM(graph, eid);
        if (VECTOR(outdegrees)[k] <= 0) {
          /* Already removed, continue */
          continue;
        }
        VECTOR(outdegrees)[k]--;
        VECTOR(outstrengths)[k] -= (weights ? VECTOR(*weights)[eid] : 1.0);
        if (VECTOR(outdegrees)[k] == 0 && VECTOR(indegrees)[k] > 0)
          IGRAPH_CHECK(igraph_dqueue_push(&sinks, k));
      }

      VECTOR(outdegrees)[v] = -1;
      VECTOR(indegrees)[v] = -1;
      nodes_left--;
    }
  }

  igraph_dqueue_destroy(&sinks);
  igraph_dqueue_destroy(&sources);
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&outstrengths);
  igraph_vector_destroy(&instrengths);
  igraph_vector_destroy(&outdegrees);
  igraph_vector_destroy(&indegrees);
  IGRAPH_FINALLY_CLEAN(7);

  /* Tidy up the ordering */
  for (i = 0; i < no_of_nodes; i++) {
    if (ordering[i] < 0)
      ordering[i] += no_of_nodes;
  }

  /* Find the feedback edges based on the ordering */
  igraph_vector_clear(result);
  j = igraph_ecount(graph);
  for (i = 0; i < j; i++) {
    if (ordering[(long)IGRAPH_FROM(graph, i)] > ordering[(long)IGRAPH_TO(graph, i)])
      IGRAPH_CHECK(igraph_vector_push_back(result, i));
  }

  igraph_free(ordering);

  IGRAPH_FINALLY_CLEAN(1);

  return IGRAPH_SUCCESS;
}

/**
 * Solves the feedback arc set problem using integer programming.
 */
int igraph_i_feedback_arc_set_ip(const igraph_t *graph, igraph_vector_t *result,
        const igraph_vector_t *weights) {
#ifndef HAVE_GLPK
  IGRAPH_ERROR("GLPK is not available", IGRAPH_UNIMPLEMENTED);    
#else
  IGRAPH_ERROR("TODO", IGRAPH_UNIMPLEMENTED);
#endif
}

