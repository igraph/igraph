/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2020  The igraph development team
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_OPERATORS_MISC_INTERNAL_H
#define IGRAPH_OPERATORS_MISC_INTERNAL_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

typedef enum {
   IGRAPH_MERGE_MODE_UNION = 1,
   IGRAPH_MERGE_MODE_INTERSECTION = 2
} igraph_i_merge_mode_t;

int igraph_i_order_edgelist_cmp(void *edges, const void *e1, const void *e2);
igraph_error_t igraph_i_merge(igraph_t *res, igraph_i_merge_mode_t mode,
                   const igraph_t *left, const igraph_t *right,
                   igraph_vector_int_t *edge_map1, igraph_vector_int_t *edge_map2);

__END_DECLS

#endif
