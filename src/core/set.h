/*
   igraph library.
   Copyright (C) 2009-2024  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CORE_SET_H
#define IGRAPH_CORE_SET_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"

IGRAPH_BEGIN_C_DECLS

/* -------------------------------------------------- */
/* Flexible set                                       */
/* -------------------------------------------------- */

/**
 * Set containing integer numbers regardless of the order
 * \ingroup types
 */

typedef struct s_set {
    igraph_int_t* stor_begin;
    igraph_int_t* stor_end;
    igraph_int_t* end;
} igraph_set_t;

#define IGRAPH_SET_NULL { 0,0,0 }
#define IGRAPH_SET_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_set_init(v, size)); \
        IGRAPH_FINALLY(igraph_set_destroy, v); } while (0)

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_set_init(igraph_set_t *set, igraph_int_t capacity);
IGRAPH_PRIVATE_EXPORT void igraph_set_destroy(igraph_set_t *set);
IGRAPH_PRIVATE_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_set_inited(igraph_set_t *set);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_set_reserve(igraph_set_t *set, igraph_int_t capacity);
IGRAPH_PRIVATE_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_set_empty(const igraph_set_t *set);
IGRAPH_PRIVATE_EXPORT void igraph_set_clear(igraph_set_t *set);
IGRAPH_PRIVATE_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_set_size(const igraph_set_t *set);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_set_add(igraph_set_t *set, igraph_int_t e);
IGRAPH_PRIVATE_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_set_contains(const igraph_set_t *set, igraph_int_t e);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_set_iterate(const igraph_set_t *set,
                                                       igraph_int_t *state,
                                                       igraph_int_t *element);

IGRAPH_END_C_DECLS

#endif
