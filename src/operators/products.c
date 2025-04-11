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
      IGRAPH_ERROR("Cartesian product of directed and undirected graphs is not supported.", IGRAPH_EINVAL);
   }

   igraph_bool_t directed = directed1 && directed2;

   igraph_integer_t vg1 = igraph_vcount(g1);
   igraph_integer_t vg2 = igraph_vcount(g2);
   igraph_integer_t vres;

   IGRAPH_SAFE_MULT(vg1, vg2, &vres);

   igraph_vector_int_t edges;
   IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);


   igraph_integer_t eg1 = igraph_ecount(g1);
   igraph_integer_t eg2 = igraph_ecount(g2);

   igraph_integer_t from, to;
   igraph_integer_t i, j;
   
   // edge (i, j) with i from g1, and j from g2
   //   will have new vertex id: i * vg2 + j

   // Edges from g1
   for (i = 0; i < eg1; ++i) {
      igraph_edge(g1, i, &from, &to);
      // for all edges (from, to) in g1, add edge from (from, j) to (to, j)
      //    for all vertex j in g2
      for (j = 0; j < vg2; ++j) {

         // SAFE MULT and SAFE ADD not needed as < vres
         igraph_integer_t v1 = from * vg2 + j;
         igraph_integer_t v2 = to * vg2 + j;

         IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v1));
         IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v2));
      }
   }

   // Edges from g2
   for (i = 0; i < eg2; ++i) {
      igraph_edge(g2, i, &from, &to);
      // for all edges (from, to) in g2, add edge from (j, from) to (j, to)
      //    for all vertex j in g1
      for (j = 0; j < vg1; ++j) {
         igraph_integer_t v1 = j * vg2 + from;
         igraph_integer_t v2 = j * vg2 + to;

         IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v1));
         IGRAPH_CHECK(igraph_vector_int_push_back(&edges, v2));
      }
   }

   IGRAPH_CHECK(igraph_create(res, &edges, vres, directed));
   igraph_vector_int_destroy(&edges);
   IGRAPH_FINALLY_CLEAN(1);

   return IGRAPH_SUCCESS;
}

igraph_error_t igraph_product(igraph_t *res,
                        const igraph_t *g1, const igraph_t *g2,
                        igraph_product_t type) {

   switch (type) {
      case IGRAPH_PRODUCT_CARTESIAN:
         IGRAPH_CHECK(cartesian_product(res, g1, g2));
         return IGRAPH_SUCCESS;
      
      default:
         IGRAPH_ERROR("Unknown graph product type.", IGRAPH_EINVAL);
   }
}
