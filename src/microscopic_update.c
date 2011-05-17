/* -*- mode: C -*-  */
/*
  Microscopic update rules for dealing with agent-level strategy revision.
  Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>

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

#include "igraph_iterators.h"
#include "igraph_interface.h"
#include "igraph_microscopic_update.h"
#include "igraph_nongraph.h"
#include "igraph_random.h"

/*
 * Internal use only.
 * Compute the cumulative proportionate values of a vector. The vector is
 * assumed to hold values associated with vertices.
 *
 * \param graph The graph object representing the game network. No error
 *        checks will be performed on this graph. You are responsible for
 *        ensuring that this is a valid graph for the particular
 *        microscopic update rule at hand.
 * \param U A vector of vertex values for which we want to compute cumulative
 *        proportionate values. The vector could be, for example, a vector of
 *        fitness for vertices of \p graph. It is assumed that each value of U
 *        is nonnegative; it is your responsibility to ensure this. Also U, or
 *        a combination of interest, is assumed to sum to a positive value;
 *        this condition will be checked.
 * \param V Pointer to an uninitialized vector. The cumulative proportionate
 *        values will be computed and stored here. No error checks will be
 *        performed on this parameter.
 * \param islocal Boolean; this flag controls which perspective to use. If
 *        true then we use the local perspective; otherwise we use the global
 *        perspective. The local perspective for a vertex v is the set of all
 *        immediate neighbours of v. In contrast, the global perspective
 *        for v is the vertex set of \p graph.
 * \param vid The vertex to use if we are considering a local perspective,
 *        i.e. if \p islocal is true. This vertex will be ignored if
 *        \p islocal is false. That is, if \p islocal is false then it is safe
 *        pass the value -1 here. On the other hand, if \p islocal is true then
 *        it is assumed that this is indeed a vertex of \p graph.
 * \param mode Defines the sort of neighbourhood to consider for \p vid. This
 *        is only relevant if we are considering the local perspective, i.e. if
 *        \p islocal is true. If we are considering the global perspective,
 *        then this parameter would be ignored. In other words, if \p islocal
 *        is false then it is safe to pass the value \p IGRAPH_ALL here. If
 *        \p graph is undirected, then we use all the immediate neighbours of
 *        \p vid. Thus if you know that \p graph is undirected, then it is
 *        safe to pass the value \p IGRAPH_ALL here. Supported values are:
 *        \clist
 *        \cli IGRAPH_OUT
 *          Use the out-neighbours of \p vid. This option is only relevant
 *          when \p graph is a digraph and we are considering the local
 *          perspective.
 *        \cli IGRAPH_IN
 *          Use the in-neighbours of \p vid. Again this option is only relevant
 *          when \p graph is a directed graph and we are considering the local
 *          perspective.
 *        \cli IGRAPH_ALL
 *          Use both the in- and out-neighbours of \p vid. This option is only
 *          relevant if \p graph is a digraph and we are considering a local
 *          perspective. Also use this value if \p graph is undirected or we
 *          are considering the global perspective.
 *        \endclist
 * \return Codes:
 *         \clist
 *         \cli IGRAPH_EINVAL
 *           This error code is returned in the following case: The vector
 *           \p U, or some combination of its values, sums to zero.
 *         \cli IGRAPH_SUCCESS
 *           This signal is returned if the cumulative proportionate values
 *           were successfully computed.
 *         \endclist
 *
 * Time complexity: O(2n) where n is the number of vertices in the
 * perspective of vid.
 */

int igraph_vcumulative_proportionate_values(const igraph_t *graph,
                                                  const igraph_vector_t *U,
                                                  igraph_vector_t *V,
                                                  igraph_bool_t islocal,
                                                  igraph_integer_t vid,
                                                  igraph_neimode_t mode) {
  igraph_integer_t v;
  igraph_real_t C;  /* cumulative probability */
  igraph_real_t P;  /* probability */
  igraph_real_t S;  /* sum of values */
  igraph_vit_t A;   /* all vertices in v's perspective */
  igraph_vs_t vs;
  long int i;

  /* Set the perspective. Let v be the vertex under consideration; it might */
  /* be that we want to update v's strategy. The local perspective for v */
  /* consists of its immediate neighbours. In contrast, the global */
  /* perspective for v are all the vertices in the given graph. Hence in the */
  /* global perspective, we will ignore the given vertex and the given */
  /* neighbourhood type, but instead consider all vertices in the given */
  /* graph. */
  if (islocal)
    IGRAPH_CHECK(igraph_vs_adj(&vs, vid, mode));
  else
    IGRAPH_CHECK(igraph_vs_all(&vs));
  IGRAPH_FINALLY(igraph_vs_destroy, &vs);

  /* Sum up all the values of vector U in the perspective for v. This */
  /* sum will be used in normalizing each value. If we are using a local */
  /* perspective, then we also need to consider the quantity of v in */
  /* computing the sum. */
  /* NOTE: Here we assume that each value to be summed is nonnegative, */
  /* and at least one of the values is nonzero. The behaviour resulting */
  /* from all values being zero would be division by zero later on when */
  /* we normalize each value. We check to see that the values sum to zero. */
  /* NOTE: In this function, the order in which we iterate through the */
  /* vertices of interest should be the same as the order in which we do so */
  /* in the caller function. If the caller function doesn't care about the */
  /* order of values in the resulting vector V, then there's no need to take */
  /* special notice of that order. But in some cases the order of values in */
  /* V is taken into account, for example, in roulette wheel selection. */
  S = 0.0;
  IGRAPH_CHECK(igraph_vit_create(graph, vs, &A));
  IGRAPH_FINALLY(igraph_vit_destroy, &A);
  while (!IGRAPH_VIT_END(A)) {
    v = (igraph_integer_t)IGRAPH_VIT_GET(A);
    S += (igraph_real_t)VECTOR(*U)[v];
    IGRAPH_VIT_NEXT(A);
  }
  if (islocal)
    S += (igraph_real_t)VECTOR(*U)[vid];
  /* avoid division by zero later on */
  if (S == (igraph_real_t)0.0) {
    igraph_vit_destroy(&A);
    igraph_vs_destroy(&vs);
    IGRAPH_FINALLY_CLEAN(2);
    IGRAPH_ERROR("Vector of values sums to zero", IGRAPH_EINVAL);
  }

  /* Get cumulative probability and relative value for each vertex in the */
  /* perspective of v. The vector V holds the cumulative proportionate */
  /* values of all vertices in v's perspective. The value V[0] is the */
  /* cumulative proportionate value of the first vertex in the vertex */
  /* iterator A. The value V[1] is the cumulative proportionate value of */
  /* the second vertex in the iterator A. And so on. If we are using the */
  /* local perspective, then we also need to consider the cumulative */
  /* proportionate value of v. In the case of the local perspective, we */
  /* don't need to compute and store v's cumulative proportionate value, */
  /* but we pretend that such value is appended to the vector V. */
  C = 0.0;
  i = 0;
  IGRAPH_VIT_RESET(A);
  IGRAPH_VECTOR_INIT_FINALLY(V, IGRAPH_VIT_SIZE(A));
  while (!IGRAPH_VIT_END(A)) {
    v = (igraph_integer_t)IGRAPH_VIT_GET(A);
    /* NOTE: Beware of division by zero here. This can happen if the vector */
    /* of values, or a combination of interest, sums to zero. */
    P = (igraph_real_t)VECTOR(*U)[v] / S;
    C += P;
    VECTOR(*V)[i] = C;
    i++;
    IGRAPH_VIT_NEXT(A);
  }

  igraph_vit_destroy(&A);
  igraph_vs_destroy(&vs);
  IGRAPH_FINALLY_CLEAN(2);

  return IGRAPH_SUCCESS;
}

/*
 * Internal use only.
 * A set of standard tests to be performed prior to strategy updates. The
 * tests contained in this function are common to many strategy revision
 * functions in this file. This function is meant to be invoked from within
 * a specific strategy update function in order to perform certain common
 * tests, including sanity checks and conditions under which no strategy
 * updates are necessary.
 *
 * \param graph The graph object representing the game network. This cannot
 *        be the empty or trivial graph, but must have at least two vertices
 *        and one edge. If \p graph has one vertex, then no strategy update
 *        would take place. Furthermore, if \p graph has at least two vertices
 *        but zero edges, then strategy update would also not take place.
 * \param vid The vertex whose strategy is to be updated. It is assumed that
 *        \p vid represents a vertex in \p graph. No checking is performed and
 *        it is your responsibility to ensure that \p vid is indeed a vertex
 *        of \p graph. If an isolated vertex is provided, i.e. the input
 *        vertex has degree 0, then no strategy update would take place and
 *        \p vid would retain its current strategy. Strategy update would also
 *        not take place if the local neighbourhood of \p vid are its
 *        in-neighbours (respectively out-neighbours), but \p vid has zero
 *        in-neighbours (respectively out-neighbours). Loops are ignored in
 *        computing the degree (in, out, all) of \p vid.
 * \param quantities A vector of quantities providing the quantity of each
 *        vertex in \p graph. Think of each entry of the vector as being
 *        generated by a function such as the fitness function for the game.
 *        So if the vector represents fitness quantities, then each vector
 *        entry is the fitness of some vertex. The length of this vector must
 *        be the same as the number of vertices in the vertex set of \p graph.
 * \param strategies A vector of the current strategies for the vertex
 *        population. Each strategy is identified with a nonnegative integer,
 *        whose interpretation depends on the payoff matrix of the game.
 *        Generally we use the strategy ID as a row or column index of the
 *        payoff matrix. The length of this vector must be the same as the
 *        number of vertices in the vertex set of \p graph.
 * \param mode Defines the sort of neighbourhood to consider for \p vid. If
 *        \p graph is undirected, then we use all the immediate neighbours of
 *        \p vid. Thus if you know that \p graph is undirected, then it is safe
 *        to pass the value \p IGRAPH_ALL here. Supported values are:
 *        \clist
 *        \cli IGRAPH_OUT
 *          Use the out-neighbours of \p vid. This option is only relevant
 *          when \p graph is a directed graph.
 *        \cli IGRAPH_IN
 *          Use the in-neighbours of \p vid. Again this option is only relevant
 *          when \p graph is a directed graph.
 *        \cli IGRAPH_ALL
 *          Use both the in- and out-neighbours of \p vid. This option is only
 *          relevant if \p graph is a digraph. Also use this value if
 *          \p graph is undirected.
 *        \endclist
 * \param updates Boolean; at the end of this test suite, this flag
 *        indicates whether to proceed with strategy revision. If true then
 *        strategy revision should proceed; otherwise there is no need to
 *        continue with revising a vertex's strategy. A caller function that
 *        invokes this function would use the value of \p updates to
 *        determine whether to proceed with strategy revision.
 * \param islocal Boolean; this flag controls which perspective to use. If
 *        true then we use the local perspective; otherwise we use the global
 *        perspective. The local perspective for \p vid is the set of all
 *        immediate neighbours of \p vid. In contrast, the global perspective
 *        for \p vid is the vertex set of \p graph.
 * \return Codes:
 *         \clist
 *         \cli IGRAPH_EINVAL
 *           This error code is returned in each of the following cases:
 *           (1) Any of the parameters \p graph, \p quantities, or
 *           \p strategies is a null pointer. (2) The vector \p quantities
 *           or \p strategies has a length different from the number of
 *           vertices in \p graph. (3) The parameter \p graph is the empty
 *           or null graph, i.e. the graph with zero vertices and edges.
 *         \cli IGRAPH_SUCCESS
 *           This signal is returned if no errors were raised. You should use
 *           the value of the boolean \p updates to decide whether to go
 *           ahead with updating a vertex's strategy.
 *         \endclist
 */

int igraph_microscopic_standard_tests(const igraph_t *graph,
                                      igraph_integer_t vid,
                                      const igraph_vector_t *quantities,
                                      const igraph_vector_t *strategies,
                                      igraph_neimode_t mode,
				      igraph_bool_t *updates, 
				      igraph_bool_t islocal) {

  igraph_integer_t nvert;
  igraph_vector_t degv;
  *updates=1;

  /* sanity checks */
  if (graph == NULL) {
    IGRAPH_ERROR("Graph is a null pointer", IGRAPH_EINVAL);
  }
  if (quantities == NULL) {
    IGRAPH_ERROR("Quantities vector is a null pointer", IGRAPH_EINVAL);
  }
  if (strategies == NULL) {
    IGRAPH_ERROR("Strategies vector is a null pointer", IGRAPH_EINVAL);
  }

  /* the empty graph */
  nvert=igraph_vcount(graph);
  if (nvert < 1) {
    IGRAPH_ERROR("Graph cannot be the empty graph", IGRAPH_EINVAL);
  }
  /* invalid vector length */
  if (nvert != (igraph_integer_t)igraph_vector_size(quantities)) {
    IGRAPH_ERROR("Size of quantities vector different from number of vertices",
                 IGRAPH_EINVAL);
  }
  if (nvert != (igraph_integer_t)igraph_vector_size(strategies)) {
    IGRAPH_ERROR("Size of strategies vector different from number of vertices",
                 IGRAPH_EINVAL);
  }

  /* Various conditions under which no strategy updates will take place. That
   * is, the vertex retains its current strategy.
   */
  /* given graph has < 2 vertices */
  if (nvert < 2) {
    *updates=0;
  }
  /* graph has >= 2 vertices, but no edges */
  if (igraph_ecount(graph) < 1) {
    *updates=0;
  }

  /* Test for vertex isolation, depending on the perspective given. For
   * undirected graphs, a given vertex v is isolated if its degree is zero.
   * If we are considering in-neighbours (respectively out-neighbours), then
   * we say that v is isolated if its in-degree (respectively out-degree) is
   * zero. In general, this vertex isolation test is only relevant if we are
   * using a local perspective, i.e. if we only consider the immediate
   * neighbours (local perspective) of v as opposed to all vertices in the
   * vertex set of the graph (global perspective).
   */
  if (islocal) {
    /* Moving on ahead with vertex isolation test, since local perspective */
    /* is requested. */
    IGRAPH_VECTOR_INIT_FINALLY(&degv, 1);
    IGRAPH_CHECK(igraph_degree(graph, &degv, igraph_vss_1(vid),
			       mode, IGRAPH_NO_LOOPS));
    if (VECTOR(degv)[0] < 1)
      *updates = 0;
    igraph_vector_destroy(&degv);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return IGRAPH_SUCCESS;
}

/**
 * \ingroup spatialgames
 * \function igraph_deterministic_optimal_imitation
 * \brief Adopt a strategy via deterministic optimal imitation.
 *
 * A simple deterministic imitation strategy where a vertex revises its
 * strategy to that which yields a local optimal. Here "local" is with
 * respect to the immediate neighbours of the vertex. The vertex retains its
 * current strategy where this strategy yields a locally optimal quantity.
 * The quantity in this case could be a measure such as fitness.
 *
 * \param graph The graph object representing the game network. This cannot
 *        be the empty or trivial graph, but must have at least two vertices
 *        and one edge. If \p graph has one vertex, then no strategy update
 *        would take place. Furthermore, if \p graph has at least two vertices
 *        but zero edges, then strategy update would also not take place.
 * \param vid The vertex whose strategy is to be updated. It is assumed that
 *        \p vid represents a vertex in \p graph. No checking is performed and
 *        it is your responsibility to ensure that \p vid is indeed a vertex
 *        of \p graph. If an isolated vertex is provided, i.e. the input
 *        vertex has degree 0, then no strategy update would take place and
 *        \p vid would retain its current strategy. Strategy update would also
 *        not take place if the local neighbourhood of \p vid are its
 *        in-neighbours (respectively out-neighbours), but \p vid has zero
 *        in-neighbours (respectively out-neighbours). Loops are ignored in
 *        computing the degree (in, out, all) of \p vid.
 * \param optimality Logical; controls the type of optimality to be used.
 *        Supported values are:
 *        \clist
 *        \cli IGRAPH_MAXIMUM
 *          Use maximum deterministic imitation, where the strategy of the
 *          vertex with maximum quantity (e.g. fitness) would be adopted. We
 *          update the strategy of \p vid to that which yields a local
 *          maximum.
 *        \cli IGRAPH_MINIMUM
 *          Use minimum deterministic imitation. That is, the strategy of the
 *          vertex with minimum quantity would be imitated. In other words,
 *          update to the strategy that yields a local minimum.
 *        \endclist
 * \param quantities A vector of quantities providing the quantity of each
 *        vertex in \p graph. Think of each entry of the vector as being
 *        generated by a function such as the fitness function for the game.
 *        So if the vector represents fitness quantities, then each vector
 *        entry is the fitness of some vertex. The length of this vector must
 *        be the same as the number of vertices in the vertex set of \p graph.
 * \param strategies A vector of the current strategies for the vertex
 *        population. The updated strategy for \p vid would be stored here.
 *        Each strategy is identified with a nonnegative integer, whose
 *        interpretation depends on the payoff matrix of the game. Generally
 *        we use the strategy ID as a row or column index of the payoff
 *        matrix. The length of this vector must be the same as the number of
 *        vertices in the vertex set of \p graph.
 * \param mode Defines the sort of neighbourhood to consider for \p vid. If
 *        \p graph is undirected, then we use all the immediate neighbours of
 *        \p vid. Thus if you know that \p graph is undirected, then it is safe
 *        to pass the value \p IGRAPH_ALL here. Supported values are:
 *        \clist
 *        \cli IGRAPH_OUT
 *          Use the out-neighbours of \p vid. This option is only relevant
 *          when \p graph is a directed graph.
 *        \cli IGRAPH_IN
 *          Use the in-neighbours of \p vid. Again this option is only relevant
 *          when \p graph is a directed graph.
 *        \cli IGRAPH_ALL
 *          Use both the in- and out-neighbours of \p vid. This option is only
 *          relevant if \p graph is a digraph. Also use this value if
 *          \p graph is undirected.
 *        \endclist
 * \return The error code \p IGRAPH_EINVAL is returned in each of the
 *         following cases: (1) Any of the parameters \p graph, \p quantities,
 *         or \p strategies is a null pointer. (2) The vector \p quantities
 *         or \p strategies has a length different from the number of vertices
 *         in \p graph. (3) The parameter \p graph is the empty or null graph,
 *         i.e. the graph with zero vertices and edges.
 *
 * Time complexity: O(2d), where d is the degree of the vertex \p vid.
 *
 * \example examples/simple/igraph_deterministic_optimal_imitation.c
 */

int igraph_deterministic_optimal_imitation(const igraph_t *graph,
                                           igraph_integer_t vid,
                                           igraph_optimal_t optimality,
                                           const igraph_vector_t *quantities,
                                           igraph_vector_t *strategies,
                                           igraph_neimode_t mode) {
  igraph_integer_t i, k, v;
  igraph_real_t q;
  igraph_vector_t adj;
  igraph_bool_t updates;

  IGRAPH_CHECK(igraph_microscopic_standard_tests(graph, vid, quantities,
                                                 strategies, mode, &updates,
                                                 /*is local?*/ 1));
  if (!updates) { return IGRAPH_SUCCESS; } /* Nothing to do */
  
  /* Choose a locally optimal strategy to imitate. This can be either maximum
   * or minimum deterministic imitation. By now we know that the given vertex v
   * has degree >= 1 and at least 1 edge. Then within its immediate
   * neighbourhood adj(v) and including v itself, there exists a vertex whose
   * strategy yields a local optimal quantity.
   */
  /* Random permutation of adj(v). This ensures that if there are multiple */
  /* candidates with an optimal strategy, then we choose one such candidate */
  /* at random. */
  IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
  IGRAPH_CHECK(igraph_neighbors(graph, &adj, vid, mode));
  IGRAPH_CHECK(igraph_fisher_yates_shuffle(&adj));
  /* maximum deterministic imitation */
  i = vid;
  q = (igraph_real_t)VECTOR(*quantities)[vid];
  if (optimality == IGRAPH_MAXIMUM) {
    for (k = 0; k < igraph_vector_size(&adj); k++) {
      v = VECTOR(adj)[k];
      if ((igraph_real_t)VECTOR(*quantities)[v] > q) {
	i = v;
	q = (igraph_real_t)VECTOR(*quantities)[v];
      }
    }
  } else { /* minimum deterministic imitation */
    for (k = 0; k < igraph_vector_size(&adj); k++) {
      v = VECTOR(adj)[k];
      if ((igraph_real_t)VECTOR(*quantities)[v] < q) {
	i = v;
	q = (igraph_real_t)VECTOR(*quantities)[v];
      }
    }
  }
  /* Now i is a vertex with a locally optimal quantity, the value of which */
  /* is q. Update the strategy of vid to that of i. */
  VECTOR(*strategies)[vid] = VECTOR(*strategies)[i];
  igraph_vector_destroy(&adj);
  IGRAPH_FINALLY_CLEAN(1);

  return IGRAPH_SUCCESS;
}

/**
 * \ingroup spatialgames
 * \function igraph_roulette_wheel_imitation
 * \brief Adopt a strategy via roulette wheel selection.
 *
 * A simple stochastic imitation strategy where a vertex revises its
 * strategy to that of a vertex u chosen proportionate to u's quantity
 * (e.g. fitness). This is a special case of stochastic imitation, where a
 * candidate is not chosen uniformly at random but proportionate to its
 * quantity.
 *
 * \param graph The graph object representing the game network. This cannot
 *        be the empty or trivial graph, but must have at least two vertices
 *        and one edge. If \p graph has one vertex, then no strategy update
 *        would take place. Furthermore, if \p graph has at least two vertices
 *        but zero edges, then strategy update would also not take place.
 * \param vid The vertex whose strategy is to be updated. It is assumed that
 *        \p vid represents a vertex in \p graph. No checking is performed and
 *        it is your responsibility to ensure that \p vid is indeed a vertex
 *        of \p graph. If an isolated vertex is provided, i.e. the input
 *        vertex has degree 0, then no strategy update would take place and
 *        \p vid would retain its current strategy. Strategy update would also
 *        not take place if the local neighbourhood of \p vid are its
 *        in-neighbours (respectively out-neighbours), but \p vid has zero
 *        in-neighbours (respectively out-neighbours). Loops are ignored in
 *        computing the degree (in, out, all) of \p vid.
 * \param islocal Boolean; this flag controls which perspective to use in
 *        computing the relative quantity. If true then we use the local
 *        perspective; otherwise we use the global perspective. The local
 *        perspective for \p vid is the set of all immediate neighbours of
 *        \p vid. In contrast, the global perspective for \p vid is the
 *        vertex set of \p graph.
 * \param quantities A vector of quantities providing the quantity of each
 *        vertex in \p graph. Think of each entry of the vector as being
 *        generated by a function such as the fitness function for the game.
 *        So if the vector represents fitness quantities, then each vector
 *        entry is the fitness of some vertex. The length of this vector must
 *        be the same as the number of vertices in the vertex set of \p graph.
 *        For the purpose of roulette wheel selection, each vector entry is
 *        assumed to be nonnegative; no checks will be performed for this. It
 *        is your responsibility to ensure that at least one entry is nonzero.
 *        Furthermore, this vector cannot be a vector of zeros; this condition
 *        will be checked.
 * \param strategies A vector of the current strategies for the vertex
 *        population. The updated strategy for \p vid would be stored here.
 *        Each strategy is identified with a nonnegative integer, whose
 *        interpretation depends on the payoff matrix of the game. Generally
 *        we use the strategy ID as a row or column index of the payoff
 *        matrix. The length of this vector must be the same as the number of
 *        vertices in the vertex set of \p graph.
 * \param mode Defines the sort of neighbourhood to consider for \p vid. This
 *        is only relevant if we are considering the local perspective, i.e. if
 *        \p islocal is true. If we are considering the global perspective,
 *        then it is safe to pass the value \p IGRAPH_ALL here. If \p graph is
 *        undirected, then we use all the immediate neighbours of \p vid. Thus
 *        if you know that \p graph is undirected, then it is safe to pass the
 *        value \p IGRAPH_ALL here. Supported values are:
 *        \clist
 *        \cli IGRAPH_OUT
 *          Use the out-neighbours of \p vid. This option is only relevant
 *          when \p graph is a digraph and we are considering the local
 *          perspective.
 *        \cli IGRAPH_IN
 *          Use the in-neighbours of \p vid. Again this option is only relevant
 *          when \p graph is a directed graph and we are considering the local
 *          perspective.
 *        \cli IGRAPH_ALL
 *          Use both the in- and out-neighbours of \p vid. This option is only
 *          relevant if \p graph is a digraph. Also use this value if
 *          \p graph is undirected or we are considering the global
 *          perspective.
 *        \endclist
 * \return The error code \p IGRAPH_EINVAL is returned in each of the following
 *         cases: (1) Any of the parameters \p graph, \p quantities, or
 *         \p strategies is a null pointer. (2) The vector \p quantities or
 *         \p strategies has a length different from the number of vertices
 *         in \p graph. (3) The parameter \p graph is the empty or null graph,
 *         i.e. the graph with zero vertices and edges. (4) The vector
 *         \p quantities sums to zero.
 *
 * Time complexity: O(n) where n is the number of vertices in the perspective
 * to consider. If we consider the global perspective, then n is the number
 * of vertices in the vertex set of \p graph. On the other hand, for the local
 * perspective n is the degree of \p vid, excluding loops.
 *
 * </para><para>
 * Reference:
 * \clist
 * \cli (Yu &amp; Gen 2010)
 *   X. Yu and M. Gen. \emb Introduction to Evolutionary Algorithms. \eme
 *   Springer, 2010, pages 18--20.
 * \endclist
 *
 * \example examples/simple/igraph_roulette_wheel_imitation.c
 */

int igraph_roulette_wheel_imitation(const igraph_t *graph,
                                    igraph_integer_t vid,
                                    igraph_bool_t islocal,
                                    const igraph_vector_t *quantities,
                                    igraph_vector_t *strategies,
                                    igraph_neimode_t mode) {
  igraph_bool_t updates;
  igraph_integer_t u;
  igraph_real_t r;    /* random number */
  igraph_vector_t V;  /* vector of cumulative proportionate quantities */
  igraph_vit_t A;     /* all vertices in v's perspective */
  igraph_vs_t vs;
  long int i;

  IGRAPH_CHECK(igraph_microscopic_standard_tests(graph, vid, quantities,
                                                 strategies, mode, &updates,
                                                 islocal));
  if (!updates)
    return IGRAPH_SUCCESS;  /* nothing further to do */

  /* set the perspective */
  if (islocal)
    IGRAPH_CHECK(igraph_vs_adj(&vs, vid, mode));
  else
    IGRAPH_CHECK(igraph_vs_all(&vs));
  IGRAPH_FINALLY(igraph_vs_destroy, &vs);
  IGRAPH_CHECK(igraph_vit_create(graph, vs, &A));
  IGRAPH_FINALLY(igraph_vit_destroy, &A);

  IGRAPH_CHECK(igraph_vcumulative_proportionate_values(graph, quantities, &V,
                                                       islocal, vid, mode));

  /* Finally, choose a vertex u to imitate. The vertex u is chosen */
  /* proportionate to its quantity. In the case of a local perspective, we */
  /* pretend that v's cumulative proportionate quantity has been appended to */
  /* the vector V. Let V be of length n so that V[n-1] is the last element */
  /* of V, and let r be a real number chosen uniformly at random from the */
  /* unit interval [0,1]. If r > V[i] for all i < n, then v defaults to */
  /* retaining its current strategy. Similarly in the case of the global */
  /* perspective, if r > V[i] for all i < n - 1 then v would adopt the */
  /* strategy of the vertex whose cumulative proportionate quantity is */
  /* V[n-1]. */
  /* NOTE: Here we assume that the order in which we iterate through the */
  /* vertices in A is the same as the order in which we do so in the */
  /* invoked function igraph_vcumulative_proportionate_values(). */
  /* Otherwise we would incorrectly associate each V[i] with a vertex in A. */
  RNG_BEGIN();
  r = RNG_UNIF01();
  RNG_END();
  i = 0;
  while (!IGRAPH_VIT_END(A)) {
    if (r <= VECTOR(V)[i]) {
      /* We have found our candidate vertex for imitation. Update strategy */
      /* of v to that of u, and exit the selection loop. */
      u = (igraph_integer_t)IGRAPH_VIT_GET(A);
      VECTOR(*strategies)[vid] = VECTOR(*strategies)[u];
      break;
    }
    i++;
    IGRAPH_VIT_NEXT(A);
  }

  /* By now, vertex v should either retain its current strategy or it has */
  /* adopted the strategy of a vertex in its perspective. Nothing else to */
  /* do, but clean up. */
  igraph_vector_destroy(&V);
  igraph_vit_destroy(&A);
  igraph_vs_destroy(&vs);
  IGRAPH_FINALLY_CLEAN(3);

  return IGRAPH_SUCCESS;
}

/**
 * \ingroup spatialgames
 * \function igraph_stochastic_imitation
 * \brief Adopt a strategy via stochastic imitation with uniform selection.
 *
 * A simple stochastic imitation strategy where a vertex revises its
 * strategy to that of a vertex chosen uniformly at random from its local
 * neighbourhood. This is called stochastic imitation via uniform selection,
 * where the strategy to imitate is chosen via some random process. For the
 * purposes of this function, we use uniform selection from a pool of
 * candidates.
 *
 * \param graph The graph object representing the game network. This cannot
 *        be the empty or trivial graph, but must have at least two vertices
 *        and one edge. If \p graph has one vertex, then no strategy update
 *        would take place. Furthermore, if \p graph has at least two vertices
 *        but zero edges, then strategy update would also not take place.
 * \param vid The vertex whose strategy is to be updated. It is assumed that
 *        \p vid represents a vertex in \p graph. No checking is performed and
 *        it is your responsibility to ensure that \p vid is indeed a vertex
 *        of \p graph. If an isolated vertex is provided, i.e. the input
 *        vertex has degree 0, then no strategy update would take place and
 *        \p vid would retain its current strategy. Strategy update would also
 *        not take place if the local neighbourhood of \p vid are its
 *        in-neighbours (respectively out-neighbours), but \p vid has zero
 *        in-neighbours (respectively out-neighbours). Loops are ignored in
 *        computing the degree (in, out, all) of \p vid.
 * \param algo This flag controls which algorithm to use in stochastic
 *        imitation. Supported values are:
 *        \clist
 *        \cli IGRAPH_IMITATE_AUGMENTED
 *          Augmented imitation. Vertex \p vid imitates the strategy of the
 *          chosen vertex u provided that doing so would increase the
 *          quantity (e.g. fitness) of \p vid. Augmented imitation can be
 *          thought of as "imitate if better".
 *        \cli IGRAPH_IMITATE_BLIND
 *          Blind imitation. Vertex \p vid blindly imitates the strategy of
 *          the chosen vertex u, regardless of whether doing so would
 *          increase or decrease the quantity of \p vid.
 *        \cli IGRAPH_IMITATE_CONTRACTED
 *          Contracted imitation. Here vertex \p vid imitates the strategy of
 *          the chosen vertex u if doing so would decrease the quantity of
 *          \p vid. Think of contracted imitation as "imitate if worse".
 *        \endclist
 * \param quantities A vector of quantities providing the quantity of each
 *        vertex in \p graph. Think of each entry of the vector as being
 *        generated by a function such as the fitness function for the game.
 *        So if the vector represents fitness quantities, then each vector
 *        entry is the fitness of some vertex. The length of this vector must
 *        be the same as the number of vertices in the vertex set of \p graph.
 * \param strategies A vector of the current strategies for the vertex
 *        population. The updated strategy for \p vid would be stored here.
 *        Each strategy is identified with a nonnegative integer, whose
 *        interpretation depends on the payoff matrix of the game. Generally
 *        we use the strategy ID as a row or column index of the payoff
 *        matrix. The length of this vector must be the same as the number of
 *        vertices in the vertex set of \p graph.
 * \param mode Defines the sort of neighbourhood to consider for \p vid. If
 *        \p graph is undirected, then we use all the immediate neighbours of
 *        \p vid. Thus if you know that \p graph is undirected, then it is safe
 *        to pass the value \p IGRAPH_ALL here. Supported values are:
 *        \clist
 *        \cli IGRAPH_OUT
 *          Use the out-neighbours of \p vid. This option is only relevant
 *          when \p graph is a directed graph.
 *        \cli IGRAPH_IN
 *          Use the in-neighbours of \p vid. Again this option is only relevant
 *          when \p graph is a directed graph.
 *        \cli IGRAPH_ALL
 *          Use both the in- and out-neighbours of \p vid. This option is only
 *          relevant if \p graph is a digraph. Also use this value if
 *          \p graph is undirected.
 *        \endclist
 * \return The error code \p IGRAPH_EINVAL is returned in each of the following
 *         cases: (1) Any of the parameters \p graph, \p quantities, or
 *         \p strategies is a null pointer. (2) The vector \p quantities or
 *         \p strategies has a length different from the number of vertices
 *         in \p graph. (3) The parameter \p graph is the empty or null graph,
 *         i.e. the graph with zero vertices and edges. (4) The parameter
 *         \p algo refers to an unsupported stochastic imitation algorithm.
 *
 * Time complexity: depends on the uniform random number generator, but should
 * usually be O(1).
 *
 * \example examples/simple/igraph_stochastic_imitation.c
 */

int igraph_stochastic_imitation(const igraph_t *graph,
                                igraph_integer_t vid,
                                igraph_imitate_algorithm_t algo,
                                const igraph_vector_t *quantities,
                                igraph_vector_t *strategies,
                                igraph_neimode_t mode) {
  igraph_bool_t updates;
  igraph_integer_t u;
  igraph_vector_t adj;
  int i;

  /* sanity checks */
  if (algo != IGRAPH_IMITATE_AUGMENTED &&
      algo != IGRAPH_IMITATE_BLIND &&
      algo != IGRAPH_IMITATE_CONTRACTED) {
    IGRAPH_ERROR("Unsupported stochastic imitation algorithm",
                 IGRAPH_EINVAL);
  }
  IGRAPH_CHECK(igraph_microscopic_standard_tests(graph, vid, quantities,
                                                 strategies, mode, &updates,
                                                 /*is local?*/ 1));
  if (!updates)
    return IGRAPH_SUCCESS;  /* nothing more to do */

  /* immediate neighbours of v */
  IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
  IGRAPH_CHECK(igraph_neighbors(graph, &adj, vid, mode));

  /* Blind imitation. Let v be the vertex whose strategy we want to revise. */
  /* Choose a vertex u uniformly at random from the immediate neighbours of */
  /* v, including v itself. Then blindly update the strategy of v to that of */
  /* u, irrespective of whether doing so would increase or decrease the */
  /* quantity (e.g. fitness) of v. Here v retains its current strategy if */
  /* the chosen vertex u is indeed v itself. */
  if (algo == IGRAPH_IMITATE_BLIND) {
    IGRAPH_CHECK(igraph_vector_push_back(&adj, vid));
    RNG_BEGIN();
    i = RNG_INTEGER(0, igraph_vector_size(&adj) - 1);
    RNG_END();
    u = VECTOR(adj)[i];
    VECTOR(*strategies)[vid] = VECTOR(*strategies)[u];
  }
  /* Augmented imitation. Let v be the vertex whose strategy we want to */
  /* revise. Let f be the quantity function for the game. Choose a vertex u */
  /* uniformly at random from the immediate neighbours of v; do not include */
  /* v. Then v imitates the strategy of u if f(u) > f(v). Otherwise v */
  /* retains its current strategy. */
  else if (algo == IGRAPH_IMITATE_AUGMENTED) {
    RNG_BEGIN();
    i = RNG_INTEGER(0, igraph_vector_size(&adj) - 1);
    RNG_END();
    u = VECTOR(adj)[i];
    if (VECTOR(*quantities)[u] > VECTOR(*quantities)[vid])
      VECTOR(*strategies)[vid] = VECTOR(*strategies)[u];
  }
  /* Contracted imitation. Let v be the vertex whose strategy we want to */
  /* update and let f be the quantity function for the game. Choose a vertex */
  /* u uniformly at random from the immediate neighbours of v, excluding v */
  /* itself. Then v imitates the strategy of u provided that f(u) < f(v). */
  /* Otherwise v retains its current strategy. */
  else if (algo == IGRAPH_IMITATE_CONTRACTED) {
    RNG_BEGIN();
    i = RNG_INTEGER(0, igraph_vector_size(&adj) - 1);
    RNG_END();
    u = VECTOR(adj)[i];
    if (VECTOR(*quantities)[u] < VECTOR(*quantities)[vid])
      VECTOR(*strategies)[vid] = VECTOR(*strategies)[u];
  }

  /* clean up */
  igraph_vector_destroy(&adj);
  IGRAPH_FINALLY_CLEAN(1);

  return IGRAPH_SUCCESS;
}
