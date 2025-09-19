/*
   igraph library.
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

#ifndef IGRAPH_MISC_GRAPHICALITY_H
#define IGRAPH_MISC_GRAPHICALITY_H

#include "igraph_decls.h"
#include "igraph_graphicality.h"

#define IGRAPH_I_MULTI_EDGES_SW 0x02 /* 010, more than one edge allowed between distinct vertices */
#define IGRAPH_I_MULTI_LOOPS_SW 0x04 /* 100, more than one self-loop allowed on the same vertex   */

IGRAPH_BEGIN_C_DECLS

igraph_error_t igraph_i_edge_type_to_loops_multiple(
    igraph_edge_type_sw_t allowed_edge_type,
    igraph_bool_t *loops, igraph_bool_t *multiple);

IGRAPH_END_C_DECLS

#endif /* IGRAPH_MISC_GRAPHICALITY_H */
