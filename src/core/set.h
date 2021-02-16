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

#ifndef IGRAPH_CORE_SET_H
#define IGRAPH_CORE_SET_H

#include "igraph_decls.h"
#include "igraph_types.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Flexible set                                       */
/* -------------------------------------------------- */

/**
 * Set containing integer numbers regardless of the order
 * \ingroup types
 */

typedef struct s_set {
    igraph_integer_t* stor_begin;
    igraph_integer_t* stor_end;
    igraph_integer_t* end;
} igraph_set_t;

#define IGRAPH_SET_NULL { 0,0,0 }
#define IGRAPH_SET_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_set_init(v, size)); \
        IGRAPH_FINALLY(igraph_set_destroy, v); } while (0)

IGRAPH_PRIVATE_EXPORT int igraph_set_init(igraph_set_t* set, long int size);
IGRAPH_PRIVATE_EXPORT void igraph_set_destroy(igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_inited(igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT int igraph_set_reserve(igraph_set_t* set, long int size);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_empty(const igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT void igraph_set_clear(igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT long int igraph_set_size(const igraph_set_t* set);
IGRAPH_PRIVATE_EXPORT int igraph_set_add(igraph_set_t* v, igraph_integer_t e);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_contains(igraph_set_t* set, igraph_integer_t e);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_iterate(igraph_set_t* set, long int* state,
                                                       igraph_integer_t* element);

__END_DECLS

#endif
