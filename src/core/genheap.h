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

#include "igraph_types.h"
#include "igraph_vector.h"

typedef struct igraph_gen2wheap_t {
    /** Maximum number of items in the heap */
    igraph_integer_t max_size;

    /** The size of an individual item */
    size_t item_size;

    /** The items themselves in the heap */
    /* TODO: currently this is always allocated to have max_size */
    void *data;

    /** qsort-style comparison function used to order items */
    int (*cmp)(const void *, const void *);

    /** An integer index associated to each item in the heap; this vector is
     * always modified in tandem with \c data . Values in this vector are
     * between 0 and size-1 */
    igraph_vector_int_t index;

    /**
     * A _reverse_ index that allows O(1) lookup of the position of a given
     * value within the \c index member. Note that it uses two special values:
     * index2[i] == 0 means that \c i is not in \c index at all, while
     * index2[i] == 1 means that \c i is in \c index but it was _deactivated_.
     * The semantics of deactivation is up to the user of the data structure
     * to decide. Other than these two special values, index2[i] == j means
     * that index[j-2] == i and data[j-2] is the corresponding item in the heap
     */
    igraph_vector_int_t index2;
} igraph_gen2wheap_t;

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_gen2wheap_init(
        igraph_gen2wheap_t *h,
        int (*cmp)(const void *, const void *),
        size_t item_size, igraph_integer_t max_size
);
IGRAPH_PRIVATE_EXPORT void igraph_gen2wheap_destroy(igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_gen2wheap_size(const igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT void igraph_gen2wheap_clear(igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_gen2wheap_empty(const igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_gen2wheap_push_with_index(igraph_gen2wheap_t *h,
                                                                      igraph_integer_t idx, const void *elem);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_gen2wheap_max_size(const igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT const void *igraph_gen2wheap_max(const igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_gen2wheap_max_index(const igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_gen2wheap_has_elem(const igraph_gen2wheap_t *h, igraph_integer_t idx);
IGRAPH_PRIVATE_EXPORT igraph_bool_t igraph_gen2wheap_has_active(const igraph_gen2wheap_t *h, igraph_integer_t idx);
IGRAPH_PRIVATE_EXPORT const void *igraph_gen2wheap_get(const igraph_gen2wheap_t *h, igraph_integer_t idx);
IGRAPH_PRIVATE_EXPORT void igraph_gen2wheap_delete_max(igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT void igraph_gen2wheap_deactivate_max(igraph_gen2wheap_t *h);
IGRAPH_PRIVATE_EXPORT void igraph_gen2wheap_modify(igraph_gen2wheap_t *h, igraph_integer_t idx, const void *elem);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_gen2wheap_check(const igraph_gen2wheap_t *h);
