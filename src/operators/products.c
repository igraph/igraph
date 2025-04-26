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

#include "igraph_constructors.h"
#include "igraph_operators.h"

#include "igraph_conversion.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"


static igraph_error_t cartesian_product(igraph_t *res,
   const igraph_t *g1, const igraph_t *g2) {
      
   igraph_bool_t directed1 = igraph_is_directed(g1);
   igraph_bool_t directed2 = igraph_is_directed(g2);

   if (directed1 != directed2) {
      IGRAPH_ERROR("Cartesian product between a directed and an undirected graph is invalid. "
         "Please ensure that both graphs have the same directionality.", IGRAPH_EINVAL);
   }
  
   igraph_bool_t directed = directed1 && directed2;

   igraph_integer_t vcountg1 = igraph_vcount(g1);
   igraph_integer_t vcountg2 = igraph_vcount(g2);
   igraph_integer_t ecountg1 = igraph_ecount(g1);
   igraph_integer_t ecountg2 = igraph_ecount(g2);
   igraph_integer_t vcount;
   igraph_integer_t ecount;
   igraph_integer_t temp;

   IGRAPH_SAFE_MULT(vcountg1, vcountg2, &vcount);
   igraph_vector_int_t edges;

   // new edge count = vcountg1*e2 + vcountg2*e1
   IGRAPH_SAFE_MULT(vcountg1, ecountg2, &ecount);
   IGRAPH_SAFE_MULT(vcountg2, ecountg1, &temp);
   IGRAPH_SAFE_ADD(ecount, temp, &ecount);
   IGRAPH_SAFE_MULT(ecount, 2, &ecount);
   IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ecount);

   igraph_integer_t from, to;
   igraph_integer_t i, j;
   
   // Vertex ((i, j)) with i from g1, and j from g2
   //   will have new vertex id: i * vcountg2 + j

   igraph_integer_t edge_index = 0;
   // Edges from g1
   for (i = 0; i < ecountg1; ++i) {
      from = IGRAPH_FROM(g1, i);
      to = IGRAPH_TO(g1, i);
      // for all edges (from, to) in g1, add edge from ((from, j)) to ((to, j))
      //    for all vertex j in g2
      for (j = 0; j < vcountg2; ++j) {
         // SAFE MULT and SAFE ADD not needed as < vcount
         VECTOR(edges)[edge_index++] = from * vcountg2 + j;   // ((from, j))
         VECTOR(edges)[edge_index++] = to * vcountg2 + j;     // ((to, j))
      }
   }

   // Edges from g2
   for (i = 0; i < ecountg2; ++i) {
      from = IGRAPH_FROM(g2, i);
      to = IGRAPH_TO(g2, i);
      // for all edges (from, to) in g2, add edge from (j, from) to (j, to)
      //    for all vertex j in g1
      for (j = 0; j < vcountg1; ++j) {
         VECTOR(edges)[edge_index++] = j * vcountg2 + from; // ((j, from))
         VECTOR(edges)[edge_index++] = j * vcountg2 + to;   // ((j, to))
      }
   }

   IGRAPH_CHECK(igraph_create(res, &edges, vcount, directed));
   igraph_vector_int_destroy(&edges);
   IGRAPH_FINALLY_CLEAN(1);

   return IGRAPH_SUCCESS;
}

static igraph_error_t tensor_product(igraph_t *res,
   const igraph_t *g1, const igraph_t *g2) {
      
   igraph_bool_t directed1 = igraph_is_directed(g1);
   igraph_bool_t directed2 = igraph_is_directed(g2);

   if (directed1 != directed2) {
      IGRAPH_ERROR("Tensor product between a directed and an undirected graph is invalid. "
         "Please ensure that both graphs have the same directionality.", IGRAPH_EINVAL);
   }

   igraph_bool_t directed = directed1 && directed2;

   igraph_integer_t vcountg1 = igraph_vcount(g1);
   igraph_integer_t vcountg2 = igraph_vcount(g2);
   igraph_integer_t ecountg1 = igraph_ecount(g1);
   igraph_integer_t ecountg2 = igraph_ecount(g2);
   igraph_integer_t vcount;
   igraph_integer_t ecount;

   IGRAPH_SAFE_MULT(vcountg1, vcountg2, &vcount);
   igraph_vector_int_t edges;

   // new edge count = 2*e1*e2 if undirected else e1*e2
   IGRAPH_SAFE_MULT(ecountg1, ecountg2, &ecount);
   if (!directed) { // directed tensor product has only e1*e2 edges, see below
      IGRAPH_SAFE_MULT(ecount, 2, &ecount);
   }
   IGRAPH_SAFE_MULT(ecount, 2, &ecount);
   IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ecount);

   igraph_integer_t from1, to1, from2, to2;
   igraph_integer_t i, j;
   
   // Vertex ((i, j)) with i from g1, and j from g2
   //   will have new vertex id: i * vcountg2 + j
   igraph_integer_t edge_index = 0;

   for (i = 0; i < ecountg1; ++i) {
      from1 = IGRAPH_FROM(g1, i);
      to1 = IGRAPH_TO(g1, i);

      for (j = 0; j < ecountg2; ++j) {
         from2 = IGRAPH_FROM(g2, j);
         to2 = IGRAPH_TO(g2, j);

         // create edge between ((from1, from2)) to ((to1, to2))
         VECTOR(edges)[edge_index++] = from1 * vcountg2 + from2; // ((from1, from2))
         VECTOR(edges)[edge_index++] = to1   * vcountg2 + to2;   // ((to1, to2))

         // this cross edge is not present in directed edge
         // as to2 is not adjacent to from2, if direction is taken in account
         if (!directed) {
            VECTOR(edges)[edge_index++] = from1 * vcountg2 + to2;   // ((from1, to2))
            VECTOR(edges)[edge_index++] = to1   * vcountg2 + from2; // ((to1, from2))
         }
      }
   }

   IGRAPH_CHECK(igraph_create(res, &edges, vcount, directed));
   igraph_vector_int_destroy(&edges);
   IGRAPH_FINALLY_CLEAN(1);

   return IGRAPH_SUCCESS;
}

/**
 * \function igraph_product
 * \brief Computes the graph product of two graphs, based on the specified product type.
 * </para><para>
 * Supported graph product types are:
 * \clist
 *    \cli IGRAPH_PRODUCT_CARTESIAN
 *       Computes the Cartesian product of two graphs \c g1 and \c g2.
 * The Cartesian product of two graphs \c g1 and \c g2 is a graph \c res such that:
 * \olist
 *    \oli The vertex set of \c res is the Cartesian product of the vertex sets of g1 and g2: V(g1) x V(g2).
 * </para><para>
 *    \oli Two vertices <code>(u, v)</code> and <code>(u1, v1)</code> are adjacent in 
 *    \c res if and only if either <code>u = u1</code> and \c v is adjacent 
 *    to \c v1 in \c g2, or <code>v = v1</code> and \c u is adjacent to 
 *    \c u1 in \c g1.
 * \endolist
 * Thus, the number of vertices in \c res is |V1| x |V2|, and the number of edges in
 * \c res is |V1| × |E2| + |V2| × |E1|.
 * </para><para>
 * Time Complexity: O(|V1| × |V2| + |V1| × |E2| + |V2| × |E1|)
 *       where |V1| and |V2| are the number of vertices, and
 *       |E1| and |E2| are the number of edges in \c g1 and \c g2 respectively.
 * 
 * \endclist
 * Both graphs must be of the same type, either directed or undirected. If a product of an undirected and a
 * directed graph is required, the undirected graph can be converted to a directed graph using \ref 
 * igraph_to_directed, or the directed graph can be converted to an undirected graph using \ref 
 * igraph_to_undirected.
 * </para><para>
 * \param res Pointer to an uninitialized graph object. The result will be stored here.
 * \param g1 The first operand graph.
 * \param g2 The second operand graph.
 * \param type The type of graph product to compute.
 * 
 * \return Error code:
 *         \c IGRAPH_EINVAL if the specified \p type is unsupported or the input graphs
 *         \p g1 and \p g2 are incompatible for the requested product.
 */

igraph_error_t igraph_product(igraph_t *res,
                        const igraph_t *g1, const igraph_t *g2,
                        igraph_product_t type) {

   switch (type) {
      case IGRAPH_PRODUCT_CARTESIAN:
         IGRAPH_CHECK(cartesian_product(res, g1, g2));
         return IGRAPH_SUCCESS;
      
      case IGRAPH_PRODUCT_TENSOR:
         IGRAPH_CHECK(tensor_product(res, g1, g2));
         return IGRAPH_SUCCESS;

      default:
         IGRAPH_ERROR("Unknown graph product type.", IGRAPH_EINVAL);
   }
}
