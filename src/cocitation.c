/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include "igraph.h"
#include "memory.h"
#include "config.h"
#include <math.h>

int igraph_cocitation_real(const igraph_t *graph, igraph_matrix_t *res, 
                           igraph_vs_t vids, igraph_neimode_t mode,
                           igraph_vector_t *weights);

/**
 * \ingroup structural
 * \function igraph_cocitation
 * \brief Cocitation coupling.
 * 
 * </para><para>
 * Two vertices are cocited if there is another vertex citing both of
 * them. \ref igraph_cocitation() simply counts how many times two vertices are
 * cocited.
 * The cocitation score for each given vertex and all other vertices
 * in the graph will be calculated.
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 * 
 * Time complexity: O(|V|d^2), |V| is
 * the number of vertices in the graph,
 * d is the (maximum) degree of 
 * the vertices in the graph.
 *
 * \sa \ref igraph_bibcoupling()
 */

int igraph_cocitation(const igraph_t *graph, igraph_matrix_t *res, 
                      const igraph_vs_t vids) {
  return igraph_cocitation_real(graph, res, vids, IGRAPH_OUT, 0);
}

/**
 * \ingroup structural
 * \function igraph_bibcoupling
 * \brief Bibliographic coupling.
 * 
 * </para><para>
 * The bibliographic coupling of two vertices is the number
 * of other vertices they both cite, \ref igraph_bibcoupling() calculates
 * this.
 * The bibliographic coupling  score for each given vertex and all
 * other vertices in the graph will be calculated.
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 * 
 * Time complexity: O(|V|d^2),
 * |V| is the number of vertices in
 * the graph, d is the (maximum)
 * degree of the vertices in the graph.
 *
 * \sa \ref igraph_cocitation()
 */

int igraph_bibcoupling(const igraph_t *graph, igraph_matrix_t *res, 
                       const igraph_vs_t vids) {
  return igraph_cocitation_real(graph, res, vids, IGRAPH_IN, 0);
}

/**
 * \ingroup structural
 * \function igraph_similarity_inverse_log_weighted
 * \brief Vertex similarity based on the inverse logarithm of vertex degrees. 
 * 
 * </para><para>
 * The inverse log-weighted similarity of two vertices is the number of
 * their common neighbors, weighted by the inverse logarithm of their degrees.
 * It is based on the assumption that two vertices should be considered
 * more similar if they share a low-degree common neighbor, since high-degree
 * common neighbors are more likely to appear even by pure chance.
 *
 * </para><para>
 * Isolated vertices will have zero similarity to any other vertex.
 * Self-similarities are not calculated.
 *
 * </para><para>
 * See the following paper for more details: Lada A. Adamic and Eytan Adar:
 * Friends and neighbors on the Web. Social Networks, 25(3):211-230, 2003.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows is the same as the
 *        number of vertex ids in \p vids, the number of
 *        columns is the number of vertices in the graph.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node. Nodes
 *          will be weighted according to their in-degree.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node. Nodes
 *          will be weighted according to their out-degree.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation. Every node is weighted according to its undirected
 *          degree.
 *        \endclist
 * \return Error code:
 *         \c IGRAPH_EINVVID: invalid vertex id.
 * 
 * Time complexity: O(|V|d^2),
 * |V| is the number of vertices in
 * the graph, d is the (maximum)
 * degree of the vertices in the graph.
 */

int igraph_similarity_inverse_log_weighted(const igraph_t *graph,
  igraph_matrix_t *res, const igraph_vs_t vids, igraph_neimode_t mode) {
  igraph_vector_t weights;
  igraph_neimode_t mode0;
  long int i, no_of_nodes;

  switch (mode) {
    case IGRAPH_OUT: mode0 = IGRAPH_IN; break;
    case IGRAPH_IN: mode0 = IGRAPH_OUT; break;
    default: mode0 = IGRAPH_ALL;
  }

  no_of_nodes = igraph_vcount(graph);

  IGRAPH_VECTOR_INIT_FINALLY(&weights, no_of_nodes);
  IGRAPH_CHECK(igraph_degree(graph, &weights, igraph_vss_all(), mode0, 1));
  for (i=0; i < no_of_nodes; i++) {
    if (VECTOR(weights)[i] > 1)
      VECTOR(weights)[i] = 1.0 / log(VECTOR(weights)[i]);
  }

  IGRAPH_CHECK(igraph_cocitation_real(graph, res, vids, mode0, &weights));
  igraph_vector_destroy(&weights);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_cocitation_real(const igraph_t *graph, igraph_matrix_t *res, 
                           igraph_vs_t vids,
                           igraph_neimode_t mode,
                           igraph_vector_t *weights) {

  long int no_of_nodes=igraph_vcount(graph);
  long int from, i, j;
  igraph_bool_t *calc;
  igraph_matrix_t tmpres=IGRAPH_MATRIX_NULL;
  igraph_vector_t neis=IGRAPH_VECTOR_NULL;
  igraph_vit_t vit;
  
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  calc=igraph_Calloc(no_of_nodes, igraph_bool_t);
  if (calc==0) {
    IGRAPH_ERROR("cannot calculate cocitation/bibcoupling", IGRAPH_ENOMEM);
  }  
  IGRAPH_FINALLY(free, calc);   /* TODO: hack */

  for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
    calc[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
  }
  
  IGRAPH_MATRIX_INIT_FINALLY(&tmpres, no_of_nodes, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_matrix_resize(res, IGRAPH_VIT_SIZE(vit), no_of_nodes));

  /* The result */
  
  for (from=0; from<no_of_nodes; from++) {
    igraph_real_t weight = 1;

    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, from, mode));
    if (weights) weight = VECTOR(*weights)[from];

    for (i=0; i < igraph_vector_size(&neis)-1; i++) {
      if (calc[ (long int)VECTOR(neis)[i] ]) {
        for (j=i+1; j<igraph_vector_size(&neis); j++) {
          MATRIX(tmpres, (long int)VECTOR(neis)[i], 
                 (long int)VECTOR(neis)[j]) += weight;
          MATRIX(tmpres, (long int)VECTOR(neis)[j], 
                 (long int)VECTOR(neis)[i]) += weight;
        }
      }
    }
  }

  /* Copy result */
  for (IGRAPH_VIT_RESET(vit), i=0; 
       !IGRAPH_VIT_END(vit); 
       IGRAPH_VIT_NEXT(vit), i++) {
    for (j=0; j<no_of_nodes; j++) {
      MATRIX(*res, i, j) = MATRIX(tmpres, (long int) IGRAPH_VIT_GET(vit), j);
    }
  }  
  
  /* Clean up */
  igraph_matrix_destroy(&tmpres);
  igraph_vector_destroy(&neis);
  igraph_Free(calc);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

int igraph_i_neisets_intersect(const igraph_vector_t *v1,
  const igraph_vector_t *v2, long int *len_union,
  long int *len_intersection) {
  /* ASSERT: v1 and v2 are sorted */
  long int i, j, i0, j0;
  i0 = igraph_vector_size(v1); j0 = igraph_vector_size(v2);
  *len_union = i0+j0; *len_intersection = 0;
  i = 0; j = 0;
  while (i < i0 && j < j0) {
    if (VECTOR(*v1)[i] == VECTOR(*v2)[j]) {
      (*len_intersection)++; (*len_union)--;
      i++; j++;
    } else if (VECTOR(*v1)[i] < VECTOR(*v2)[j]) i++;
    else j++;
  }
  return 0;
}

/**
 * \ingroup structural
 * \function igraph_similarity_jaccard
 * \brief Jaccard similarity coefficient.
 *
 * </para><para>
 * The Jaccard similarity coefficient of two vertices is the number of common
 * neighbors divided by the number of vertices that are neighbors of at
 * least one of the two vertices being considered. This function calculates
 * the pairwise Jaccard similarities for some (or all) of the vertices.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows and columns is the same
 *        as the number of vertex ids in \p vids.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves in the neighbor
 *        sets.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 * 
 * Time complexity: O(|V|^2 d),
 * |V| is the number of vertices in the vertex iterator given, d is the
 * (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_dice(), a measure very similar to the Jaccard
 *   coefficient
 */
int igraph_similarity_jaccard(const igraph_t *graph, igraph_matrix_t *res,
    const igraph_vs_t vids, igraph_neimode_t mode, igraph_bool_t loops) {
  igraph_lazy_adjlist_t al;
  igraph_vit_t vit, vit2;
  long int i, j, k;
  long int len_union, len_intersection;
  igraph_vector_t *v1, *v2;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit2));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit2);

  IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, mode, IGRAPH_SIMPLIFY));
  IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

  IGRAPH_CHECK(igraph_matrix_resize(res, IGRAPH_VIT_SIZE(vit), IGRAPH_VIT_SIZE(vit)));

  if (loops) {
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
      i=IGRAPH_VIT_GET(vit);
      v1=igraph_lazy_adjlist_get(&al, i);
      if (!igraph_vector_binsearch(v1, i, &k)) igraph_vector_insert(v1, k, i);
    }
  }

  for (IGRAPH_VIT_RESET(vit), i=0;
    !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    MATRIX(*res, i, i) = 1.0;
    for (IGRAPH_VIT_RESET(vit2), j=0;
      !IGRAPH_VIT_END(vit2); IGRAPH_VIT_NEXT(vit2), j++) {
      if (j <= i) continue;
      v1=igraph_lazy_adjlist_get(&al, IGRAPH_VIT_GET(vit));
      v2=igraph_lazy_adjlist_get(&al, IGRAPH_VIT_GET(vit2));
      igraph_i_neisets_intersect(v1, v2, &len_union, &len_intersection);
      if (len_union > 0)
        MATRIX(*res, i, j) = ((igraph_real_t)len_intersection)/len_union;
      else
        MATRIX(*res, i, j) = 0.0;
      MATRIX(*res, j, i) = MATRIX(*res, i, j);
    }
  }

  igraph_lazy_adjlist_destroy(&al);
  igraph_vit_destroy(&vit);
  igraph_vit_destroy(&vit2);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}


/**
 * \ingroup structural
 * \function igraph_similarity_dice
 * \brief Dice similarity coefficient.
 *
 * </para><para>
 * The Dice similarity coefficient of two vertices is twice the number of common
 * neighbors divided by the sum of the degrees of the vertices. This function
 * calculates the pairwise Dice similarities for some (or all) of the vertices.
 *
 * \param graph The graph object to analyze
 * \param res Pointer to a matrix, the result of the calculation will
 *        be stored here. The number of its rows and columns is the same
 *        as the number of vertex ids in \p vids.
 * \param vids The vertex ids of the vertices for which the
 *        calculation will be done.
 * \param mode The type of neighbors to be used for the calculation in
 *        directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing edges will be considered for each node.
 *        \cli IGRAPH_IN
 *          the incoming edges will be considered for each node.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an undirected one for the
 *          computation.
 *        \endclist
 * \param loops Whether to include the vertices themselves as their own
 *        neighbors.
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           invalid vertex id passed.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 * 
 * Time complexity: O(|V|^2 d),
 * |V| is the number of vertices in the vertex iterator given, d is the
 * (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_jaccard(), a measure very similar to the Dice
 *   coefficient
 */
int igraph_similarity_dice(const igraph_t *graph, igraph_matrix_t *res,
    const igraph_vs_t vids, igraph_neimode_t mode, igraph_bool_t loops) {
  igraph_lazy_adjlist_t al;
  igraph_vit_t vit, vit2;
  long int i, j, k;
  long int len_union, len_intersection;
  igraph_vector_t *v1, *v2;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit2));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit2);

  IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, mode, IGRAPH_SIMPLIFY));
  IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

  IGRAPH_CHECK(igraph_matrix_resize(res, IGRAPH_VIT_SIZE(vit), IGRAPH_VIT_SIZE(vit)));

  if (loops) {
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
      i=IGRAPH_VIT_GET(vit);
      v1=igraph_lazy_adjlist_get(&al, i);
      if (!igraph_vector_binsearch(v1, i, &k)) igraph_vector_insert(v1, k, i);
    }
  }

  for (IGRAPH_VIT_RESET(vit), i=0;
    !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
    MATRIX(*res, i, i) = 1.0;
    for (IGRAPH_VIT_RESET(vit2), j=0;
      !IGRAPH_VIT_END(vit2); IGRAPH_VIT_NEXT(vit2), j++) {
      if (j <= i) continue;
      v1=igraph_lazy_adjlist_get(&al, IGRAPH_VIT_GET(vit));
      v2=igraph_lazy_adjlist_get(&al, IGRAPH_VIT_GET(vit2));
      igraph_i_neisets_intersect(v1, v2, &len_union, &len_intersection);
      len_union += len_intersection;
      if (len_union > 0)
        MATRIX(*res, i, j) = (2.0*len_intersection)/len_union;
      else
        MATRIX(*res, i, j) = 0.0;
      MATRIX(*res, j, i) = MATRIX(*res, i, j);
    }
  }

  igraph_lazy_adjlist_destroy(&al);
  igraph_vit_destroy(&vit);
  igraph_vit_destroy(&vit2);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

