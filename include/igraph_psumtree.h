/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_PSUMTREE_H
#define IGRAPH_PSUMTREE_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

IGRAPH_BEGIN_C_DECLS

typedef struct {
    igraph_vector_t v;
    igraph_int_t size;
    igraph_int_t offset;
} igraph_psumtree_t;

IGRAPH_EXPORT igraph_error_t igraph_psumtree_init(igraph_psumtree_t *t, igraph_int_t size);
IGRAPH_EXPORT void igraph_psumtree_reset(igraph_psumtree_t *t);
IGRAPH_EXPORT void igraph_psumtree_destroy(igraph_psumtree_t *t);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t igraph_psumtree_get(const igraph_psumtree_t *t, igraph_int_t idx);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_psumtree_size(const igraph_psumtree_t *t);
IGRAPH_EXPORT igraph_error_t igraph_psumtree_search(const igraph_psumtree_t *t, igraph_int_t *idx,
                                         igraph_real_t elem);
IGRAPH_EXPORT igraph_error_t igraph_psumtree_update(igraph_psumtree_t *t, igraph_int_t idx,
                                         igraph_real_t new_value);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t igraph_psumtree_sum(const igraph_psumtree_t *t);

IGRAPH_END_C_DECLS

#endif
