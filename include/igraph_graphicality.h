/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2009-2020  Gabor Csardi <csardi.gabor@gmail.com>

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
#include "igraph_datatype.h"

__BEGIN_DECLS

typedef unsigned char igraph_edge_type_sw_t;

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

IGRAPH_EXPORT int igraph_is_graphical(const igraph_vector_t *out_degrees,
                                      const igraph_vector_t *in_degrees,
                                      const igraph_edge_type_sw_t allowed_edge_types,
                                      igraph_bool_t *res);

IGRAPH_EXPORT int igraph_is_bigraphical(const igraph_vector_t *degrees1,
                                        const igraph_vector_t *degrees2,
                                        const igraph_edge_type_sw_t allowed_edge_types,
                                        igraph_bool_t *res);


/* Legacy functions (deprecated): */

IGRAPH_EXPORT IGRAPH_DEPRECATED int igraph_is_degree_sequence(const igraph_vector_t *out_degrees,
                                                              const igraph_vector_t *in_degrees,
                                                              igraph_bool_t *res);

IGRAPH_EXPORT IGRAPH_DEPRECATED int igraph_is_graphical_degree_sequence(const igraph_vector_t *out_degrees,
                                                                        const igraph_vector_t *in_degrees,
                                                                        igraph_bool_t *res);

__END_DECLS

#endif // IGRAPH_GRAPHICALITY_H
