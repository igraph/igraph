/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021 The igraph development team

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

#ifndef IGRAPH_FLOW_INTERNAL_H
#define IGRAPH_FLOW_INTERNAL_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_types.h"

#include "core/estack.h"
#include "core/marked_queue.h"

__BEGIN_DECLS

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_i_all_st_cuts_pivot(
   const igraph_t *graph, const igraph_marked_queue_int_t *S,
   const igraph_estack_t *T, igraph_integer_t source, igraph_integer_t target,
   igraph_integer_t *v, igraph_vector_int_t *Isv, void *arg);

igraph_error_t igraph_i_split_vertices(const igraph_t* graph, igraph_t* result);

__END_DECLS

#endif
