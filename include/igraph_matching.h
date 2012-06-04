/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2012  Tamas Nepusz <ntamas@gmail.com>
   
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

#ifndef IGRAPH_MATCHING_H
#define IGRAPH_MATCHING_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Matchings in graphs                                */
/* -------------------------------------------------- */

int igraph_is_matching(const igraph_t* graph,
    const igraph_vector_bool_t* types, const igraph_vector_long_t* matching,
    igraph_bool_t* result);
int igraph_is_maximal_matching(const igraph_t* graph,
    const igraph_vector_bool_t* types, const igraph_vector_long_t* matching,
    igraph_bool_t* result);

int igraph_maximum_bipartite_matching(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights, igraph_real_t eps);

int igraph_maximum_matching(const igraph_t* graph, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights);

__END_DECLS

#endif
