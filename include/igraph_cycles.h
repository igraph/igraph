/*
   IGraph library.
   Copyright (C) 2022-2024  The igraph development team

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

#ifndef IGRAPH_CYCLES_H
#define IGRAPH_CYCLES_H

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"
#include "igraph_vector_list.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_fundamental_cycles(
        const igraph_t *graph,
        igraph_vector_int_list_t *result,
        igraph_integer_t start_vid,
        igraph_integer_t bfs_cutoff,
        const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_minimum_cycle_basis(
        const igraph_t *graph,
        igraph_vector_int_list_t *result,
        igraph_integer_t bfs_cutoff,
        igraph_bool_t complete,
        igraph_bool_t use_cycle_order,
        const igraph_vector_t *weights);

IGRAPH_EXPORT igraph_error_t igraph_find_cycle(
        const igraph_t *graph,
        igraph_vector_int_t *vertices,
        igraph_vector_int_t *edges,
        igraph_neimode_t mode);

__END_DECLS

#endif
