/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CACHING_H
#define IGRAPH_CACHING_H

#include "igraph_datatype.h"
#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_types.h"

#include "internal/hacks.h"

#include <string.h> /* memset */

__BEGIN_DECLS

struct igraph_i_property_cache_t {
    igraph_bool_t value[IGRAPH_PROP_I_SIZE];

    /** Bit field that stores which of the properties are cached at the moment */
    uint32_t known;
};

igraph_error_t igraph_i_property_cache_init(igraph_i_property_cache_t *cache);
igraph_error_t igraph_i_property_cache_copy(
        igraph_i_property_cache_t *cache,
        const igraph_i_property_cache_t *other_cache);
void igraph_i_property_cache_destroy(igraph_i_property_cache_t *cache);

void igraph_i_property_cache_invalidate_conditionally(
    const igraph_t *graph, uint32_t keep_always, uint32_t keep_when_false, uint32_t keep_when_true
);

__END_DECLS

#endif /* IGRAPH_CACHING_H */
