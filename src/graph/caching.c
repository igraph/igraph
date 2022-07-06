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
 * \brief Initializes a property cache, ensuring that all values are unknown.
 */
igraph_error_t igraph_i_property_cache_init(igraph_i_property_cache_t *cache) {
    memset(cache->value, 0, sizeof(cache->value) / sizeof(cache->value[0]));
    memset(cache->known, 0, sizeof(cache->known) / sizeof(cache->known[0]));
    return IGRAPH_SUCCESS;
}

/**
 * \brief Copies a property cache.
 */
igraph_error_t igraph_i_property_cache_copy(
        igraph_i_property_cache_t *cache,
        const igraph_i_property_cache_t *other_cache) {
    *cache = *other_cache;
    return IGRAPH_SUCCESS;
}

/**
 * \brief Destroys a property cache.
 */
void igraph_i_property_cache_destroy(igraph_i_property_cache_t *cache) {
    IGRAPH_UNUSED(cache);
    /* Nothing to do */
}

/***** Developer fuctions, exposed *****/

/**
 * \brief Returns the value of a cached boolean property.
 *
 * This function provides valid results only when the property is already
 * cached. Use \ref igraph_i_property_cache_has() to retrieve whether the
 * property is cached.
 *
 * \param graph  the graph whose cache is to be checked
 * \param prop   the property to retrieve from the cache
 * \return the cached value of the property if the value is in the cache, or
 *         an undefined value otherwise
 */
igraph_bool_t igraph_i_property_cache_get_bool(const igraph_t *graph, igraph_cached_property_t prop) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_SIZE);
    assert(graph->cache != NULL);
    return graph->cache->value[prop];
}

/**
 * \brief Returns whether the cache contains a value for the given cached property.
 *
 * \param graph  the graph whose cache is to be checked
 * \param prop   the property to check in the cache
 */
igraph_bool_t igraph_i_property_cache_has(const igraph_t *graph, igraph_cached_property_t prop) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_SIZE);
    assert(graph->cache != NULL);
    return graph->cache->known[prop];
}

/**
 * \brief Stores a property value in the cache.
 *
 * \param graph  the graph whose cache is to be modified
 * \param prop   the property to update in the cache
 * \param value  the value of the property to add to the cache
 */
void igraph_i_property_cache_set_bool(const igraph_t *graph, igraph_cached_property_t prop, igraph_bool_t value) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_SIZE);
    assert(graph->cache != NULL);
    /* Even though graph is const, updating the cache is not considered modification.
     * Functions that merely compute graph properties, and thus leave the graph structure
     * intact, will often update the cache. */
    graph->cache->value[prop] = value;
    graph->cache->known[prop] = 1;
}

/**
 * \brief Invalidates the cached value of a property in a graph.
 *
 * \param graph  the graph whose cache is to be modified
 * \param prop   the property to invalidate in the cache
 */
void igraph_i_property_cache_invalidate(const igraph_t *graph, igraph_cached_property_t prop) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_SIZE);
    assert(graph->cache != NULL);
    graph->cache->known[prop] = 0;
}

/**
 * \brief Invalidates all cached properties of the graph.
 *
 * This function is typically called after the graph is modified.
 *
 * \param graph  the graph whose cache is to be invalidated
 */
void igraph_i_property_cache_invalidate_all(const igraph_t *graph) {
    for (igraph_cached_property_t prop = 0; prop < IGRAPH_PROP_I_SIZE; ++prop) {
        igraph_i_property_cache_invalidate(graph, prop);
    }
}
