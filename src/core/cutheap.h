/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2020  The igraph development team

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

#ifndef IGRAPH_CORE_CUTHEAP_H
#define IGRAPH_CORE_CUTHEAP_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* Special maximum heap, needed for the minimum cut algorithm */

typedef struct igraph_i_cutheap_t {
    igraph_vector_t heap;
    igraph_vector_t index;
    igraph_vector_t hptr;
    long int dnodes;
} igraph_i_cutheap_t;

IGRAPH_PRIVATE_EXPORT int igraph_i_cutheap_init(igraph_i_cutheap_t *ch, igraph_integer_t nodes);
IGRAPH_PRIVATE_EXPORT void igraph_i_cutheap_destroy(igraph_i_cutheap_t *ch);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_i_cutheap_empty(igraph_i_cutheap_t *ch);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_i_cutheap_active_size(igraph_i_cutheap_t *ch);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_i_cutheap_size(igraph_i_cutheap_t *ch);
IGRAPH_PRIVATE_EXPORT igraph_real_t igraph_i_cutheap_maxvalue(igraph_i_cutheap_t *ch);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_i_cutheap_popmax(igraph_i_cutheap_t *ch);
IGRAPH_PRIVATE_EXPORT int igraph_i_cutheap_update(igraph_i_cutheap_t *ch, igraph_integer_t index,
                                                  igraph_real_t add);
IGRAPH_PRIVATE_EXPORT int igraph_i_cutheap_reset_undefine(igraph_i_cutheap_t *ch, long int vertex);

__END_DECLS

#endif
