/*
   igraph library.
   Copyright (C) 2020-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_GRAPHICALITY_H
#define IGRAPH_GRAPHICALITY_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

IGRAPH_BEGIN_C_DECLS

/**
 * \typedef igraph_edge_type_sw_t
 * \brief What types of non-simple edges to allow?
 *
 * This type is used with multiple functions to specify what types of non-simple
 * edges to allow, create or consider a graph. The constants below are treated
 * as "switches" that can be turned on individually and combined using the
 * bitwise-or operator. For example,
 * <code>IGRAPH_LOOPS_SW</code>
 * allows only self-loops but not multi-edges, while
 * <code>IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW</code>
 * allows both.
 *
 * \enumval IGRAPH_SIMPLE_SW A shorthand for simple graphs only, which is the default
 *    assumption.
 * \enumval IGRAPH_LOOPS_SW Allow or consider self-loops.
 * \enumval IGRAPH_MULTI_SW Allow or consider multi-edges.
 */
typedef unsigned int igraph_edge_type_sw_t;

/*
 * bit 0: self-loops alowed?
 * bit 1: more than one edge allowed between distinct vertices?
 * bit 2: more than one self-loop allowed (assuming bit 0 is set)?
 */
enum {
  IGRAPH_SIMPLE_SW = 0x00, /* 000 */
  IGRAPH_LOOPS_SW  = 0x01, /* 001 */
  IGRAPH_MULTI_SW  = 0x06  /* 110 */
};

IGRAPH_EXPORT igraph_error_t igraph_is_graphical(const igraph_vector_int_t *out_degrees,
                                      const igraph_vector_int_t *in_degrees,
                                      igraph_edge_type_sw_t allowed_edge_types,
                                      igraph_bool_t *res);

IGRAPH_EXPORT igraph_error_t igraph_is_bigraphical(const igraph_vector_int_t *degrees1,
                                        const igraph_vector_int_t *degrees2,
                                        igraph_edge_type_sw_t allowed_edge_types,
                                        igraph_bool_t *res);

IGRAPH_END_C_DECLS

#endif // IGRAPH_GRAPHICALITY_H
