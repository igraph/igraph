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

#include "igraph_interface.h"

#include "graph/caching.h"

#include <assert.h>

/****** Strictly internal functions ******/

/**
 * Initialize property cache, ensuring that all values are unknown.
 */
igraph_error_t igraph_i_property_cache_init(igraph_i_property_cache_t *cache) {
    memset(cache->value, 0, sizeof(cache->value) / sizeof(cache->value[0]));
    memset(cache->known, 0, sizeof(cache->known) / sizeof(cache->known[0]));
    return IGRAPH_SUCCESS;
}

/**
 * Copy property cache.
 */
igraph_error_t igraph_i_property_cache_copy(
        igraph_i_property_cache_t *cache,
        const igraph_i_property_cache_t *other_cache) {
    *cache = *other_cache;
    return IGRAPH_SUCCESS;
}

/**
 * Destroy property cache.
 */
void igraph_i_property_cache_destroy(igraph_i_property_cache_t *cache) {
    IGRAPH_UNUSED(cache);
    /* Nothing to do */
}

/***** Developer fuctions, exposed *****/

/**
 * Get value for given property. Valid only when cache_known() is true.
 */
igraph_bool_t igraph_i_property_cache_value(const igraph_t *graph, igraph_property_t prop) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_LAST);
    assert(graph->cache != NULL);
    return graph->cache->value[prop];
}

/**
 * Is the value know for this property?
 */
igraph_bool_t igraph_i_property_cache_known(const igraph_t *graph, igraph_property_t prop) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_LAST);
    assert(graph->cache != NULL);
    return graph->cache->known[prop];
}

/**
 * Store a property value in the cache.
 */
void igraph_i_property_cache_set(const igraph_t *graph, igraph_property_t prop, igraph_bool_t value) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_LAST);
    assert(graph->cache != NULL);
    /* Even though graph is const, updating the cache is not considered modification.
     * Functions that merely compute graph properties, and thus leave the graph structure
     * intact, will often update the cache. */
    graph->cache->value[prop] = value;
    graph->cache->known[prop] = 1;
}

/**
 * Clear the cache of a given property in this graph.
 */
void igraph_i_property_cache_clear(const igraph_t *graph, igraph_property_t prop) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_LAST);
    assert(graph->cache != NULL);
    graph->cache->known[prop] = 0;
}

/**
 * Clear all property caches in this graph.
 */
void igraph_i_property_cache_clear_all(const igraph_t *graph) {
    for (igraph_property_t prop = 0; prop < IGRAPH_PROP_I_LAST; ++prop) {
        igraph_i_property_cache_clear(graph, prop);
    }
}
