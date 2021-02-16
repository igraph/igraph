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
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

#define IGRAPH_MERGE_MODE_UNION        1
#define IGRAPH_MERGE_MODE_INTERSECTION 2

int igraph_i_order_edgelist_cmp(void *edges, const void *e1, const void *e2);
int igraph_i_merge(igraph_t *res, int mode,
                   const igraph_t *left, const igraph_t *right,
                   igraph_vector_t *edge_map1, igraph_vector_t *edge_map2);
void igraph_i_union_intersection_destroy_vectors(igraph_vector_ptr_t *v);
void igraph_i_union_intersection_destroy_vector_longs(igraph_vector_ptr_t *v);

__END_DECLS

#endif
