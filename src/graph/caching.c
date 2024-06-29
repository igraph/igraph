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
    IGRAPH_STATIC_ASSERT(IGRAPH_PROP_I_SIZE <= 32);

    memset(cache->value, 0, sizeof(cache->value));
    cache->known = 0;
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
    return graph->cache->known & (1 << prop);
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
    graph->cache->known |= (1 << prop);
}

/**
 * \brief Stores a property value in the cache.
 *
 * This function asserts that if the value of \p prop was already known,
 * then \p value is consistent with the previously stored value.
 * If this is not the case, a fatal error is triggered, with the reasoning
 * that the cache must have become invalid/inconsistent due to a bug.
 *
 * Therefore, this function cannot be used to change an already stored
 * property to a different value. If this is your intention, invalidate
 * the cache explicitly first.
 *
 * \param graph  the graph whose cache is to be modified
 * \param prop   the property to update in the cache
 * \param value  the value of the property to add to the cache
 */
void igraph_i_property_cache_set_bool_checked(const igraph_t *graph, igraph_cached_property_t prop, igraph_bool_t value) {
    IGRAPH_ASSERT(prop >= 0 && prop < IGRAPH_PROP_I_SIZE);
    assert(graph->cache != NULL);
    /* Even though graph is const, updating the cache is not considered modification.
     * Functions that merely compute graph properties, and thus leave the graph structure
     * intact, will often update the cache. */
    if (graph->cache->known & (1 << prop)) {
        IGRAPH_ASSERT(graph->cache->value[prop] == value);
    } else {
        igraph_i_property_cache_set_bool(graph, prop, value);
    }
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
    graph->cache->known &= ~(1 << prop);
}

/**
 * \brief Invalidates all cached properties of the graph.
 *
 * This function is typically called after the graph is modified.
 *
 * \param graph  the graph whose cache is to be invalidated
 */
void igraph_i_property_cache_invalidate_all(const igraph_t *graph) {
    assert(graph->cache != NULL);
    graph->cache->known = 0;
}

/**
 * \brief Invalidates all but a few cached properties of the graph, subject to specific conditions.
 *
 * This function is typically called after the graph is modified if we know that
 * the modification does not affect certain cached properties in certain cases.
 * For instance, adding more vertices does not make a connected graph disconnected,
 * so we can keep the cached properties related to graph connectivity if they
 * were already cached as true, but we need to invalidate them if they were
 * cached as false.
 *
 * </para><para>
 * Use <code>1 << IGRAPH_PROP_SOMETHING</code> to encode an individual property
 * in the bits of the bitmask used in the arguments of this function.
 *
 * \param graph       the graph whose cache is to be invalidated
 * \param keep_always bitmask where the i-th bit corresponds to cached property \em i
 *        and it should be set to 1 if the property should be \em kept ,
 *        irrespectively of its current cached value.
 */
void igraph_i_property_cache_invalidate_conditionally(
    const igraph_t *graph, uint32_t keep_always, uint32_t keep_when_false,
    uint32_t keep_when_true
) {
    uint32_t invalidate = ~keep_always;
    uint32_t mask;
    uint32_t maybe_keep;
    igraph_bool_t cached_value;

    assert(graph->cache != NULL);

    /* The bits of maybe_keep are set to 1 for those properties that are:
     *
     * - currently cached
     * - should _probably_ be invalidated
     * - _but_ the current cached value of the property may change the decision
     */
    maybe_keep = graph->cache->known & invalidate & (keep_when_false | keep_when_true);

    if (maybe_keep) {
        for (igraph_cached_property_t prop = (igraph_cached_property_t ) 0; prop < IGRAPH_PROP_I_SIZE; ++prop) {
            mask = 1 << prop;
            if (maybe_keep & mask) {
                /* if we get here, we know that the property is cached; we have
                 * masked maybe_keep with graph->cache->known */
                cached_value = igraph_i_property_cache_get_bool(graph, prop);
                if (
                    ((keep_when_false & mask) && !cached_value) ||
                    ((keep_when_true & mask) && cached_value)
                ) {
                    invalidate &= ~mask;
                }
            }
        }
    }

    graph->cache->known &= ~invalidate;
}
