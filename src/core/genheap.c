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

#include "igraph_memory.h"

#include "core/genheap.h"

#include <string.h> /* memcpy */

#if defined(_MSC_VER) && _MSC_VER < 1927
    /* MSVC does not understand restrict before version 19.27 */
    #define restrict __restrict
#endif

#define PARENT(x)     (((x)+1)/2-1)
#define LEFTCHILD(x)  (((x)+1)*2-1)
#define RIGHTCHILD(x) (((x)+1)*2)

#define ELEM(h, x) ((char *) h->data + x * h->item_size)

/* This is a smart indexed heap that can hold elements of arbitrary type.
 * In addition to the "normal" indexed heap, it allows to access every element
 * through its index in O(1) time. In other words, for this heap the indexing
 * operation is O(1), the normal heap does this in O(n) time.
 */

static void swapfunc(char * restrict a, char * restrict b, size_t es) {
    char t;

    do {
        t = *a;
        *a++ = *b;
        *b++ = t;
    } while (--es > 0);
}

static void igraph_i_gen2wheap_switch(igraph_gen2wheap_t *h,
                                      igraph_integer_t e1, igraph_integer_t e2) {
    if (e1 != e2) {
        igraph_integer_t tmp1, tmp2;

        swapfunc(ELEM(h, e1), ELEM(h, e2), h->item_size);

        tmp1 = VECTOR(h->index)[e1];
        tmp2 = VECTOR(h->index)[e2];

        VECTOR(h->index2)[tmp1] = e2 + 2;
        VECTOR(h->index2)[tmp2] = e1 + 2;

        VECTOR(h->index)[e1] = tmp2;
        VECTOR(h->index)[e2] = tmp1;
    }
}

static void igraph_i_gen2wheap_shift_up(igraph_gen2wheap_t *h,
                                        igraph_integer_t elem) {
    if (elem == 0 || h->cmp(ELEM(h, elem), ELEM(h, PARENT(elem))) < 0) {
        /* at the top */
    } else {
        igraph_i_gen2wheap_switch(h, elem, PARENT(elem));
        igraph_i_gen2wheap_shift_up(h, PARENT(elem));
    }
}

static void igraph_i_gen2wheap_sink(igraph_gen2wheap_t *h,
                                    igraph_integer_t head) {
    igraph_integer_t size = igraph_gen2wheap_size(h);
    if (LEFTCHILD(head) >= size) {
        /* no subtrees */
    } else if (RIGHTCHILD(head) == size ||
               h->cmp(ELEM(h, LEFTCHILD(head)), ELEM(h, RIGHTCHILD(head))) >= 0) {
        /* sink to the left if needed */
        if (h->cmp(ELEM(h, head), ELEM(h, LEFTCHILD(head))) < 0) {
            igraph_i_gen2wheap_switch(h, head, LEFTCHILD(head));
            igraph_i_gen2wheap_sink(h, LEFTCHILD(head));
        }
    } else {
        /* sink to the right */
        if (h->cmp(ELEM(h, head), ELEM(h, RIGHTCHILD(head))) < 0) {
            igraph_i_gen2wheap_switch(h, head, RIGHTCHILD(head));
            igraph_i_gen2wheap_sink(h, RIGHTCHILD(head));
        }
    }
}

/* ------------------ */
/* These are public   */
/* ------------------ */

/**
 * Initializes a new two-way heap. The max_size parameter defines the maximum
 * number of items that the heap can hold.
 */
igraph_error_t igraph_gen2wheap_init(
        igraph_gen2wheap_t *h,
        int (*cmp)(const void *, const void *),
        size_t item_size, igraph_integer_t max_size
) {
    /* TODO: Currently, storage is allocated for the maximum number of elements
     * right from the start. This is sufficient for the only use case as of this
     * writing, the D-SATUR graph colouring algorithm, but it may not be efficcient
     * for other use cases. Consider improving this in the future.
     */
    h->max_size = max_size;
    /* We start with the biggest */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&h->index2, max_size);
    h->cmp = cmp;
    h->item_size = item_size;
    h->data = igraph_calloc(max_size, item_size);
    IGRAPH_CHECK_OOM(h->data, "Cannot initialize generic heap.");
    IGRAPH_FINALLY(igraph_free, h->data);
    IGRAPH_CHECK(igraph_vector_int_init(&h->index, 0));

    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

/**
 * Destroys a two-way heap.
 */
void igraph_gen2wheap_destroy(igraph_gen2wheap_t *h) {
    IGRAPH_FREE(h->data);
    igraph_vector_int_destroy(&h->index);
    igraph_vector_int_destroy(&h->index2);
}

/**
 * Returns the current number of elements in the two-way heap.
 */
igraph_integer_t igraph_gen2wheap_size(const igraph_gen2wheap_t *h) {
    return igraph_vector_int_size(&h->index);
}

/**
 * Clears a two-way heap, i.e. removes all the elements from the heap.
 */
void igraph_gen2wheap_clear(igraph_gen2wheap_t *h) {
    igraph_vector_int_clear(&h->index);
    igraph_vector_int_null(&h->index2);
}

/**
 * Returns whether the two-way heap is empty.
 */
igraph_bool_t igraph_gen2wheap_empty(const igraph_gen2wheap_t *h) {
    return igraph_vector_int_empty(&h->index);
}

/**
 * Pushes a new element into the two-way heap, with the given associated index.
 * The index must be between 0 and the size of the heap minus 1, inclusive. It
 * is assumed (and not checked) that the heap does not have another item with
 * the same index.
 */
igraph_error_t igraph_gen2wheap_push_with_index(igraph_gen2wheap_t *h,
                                                igraph_integer_t idx, const void *elem) {

    igraph_integer_t size = igraph_vector_int_size(&h->index);

    if (size > IGRAPH_INTEGER_MAX - 2) {
        /* to allow size+2 below */
        IGRAPH_ERROR("Cannot push to gen2wheap, already at maximum size.", IGRAPH_EOVERFLOW);
    }

    memcpy(ELEM(h, size), elem, h->item_size);
    IGRAPH_CHECK(igraph_vector_int_push_back(&h->index, idx));
    VECTOR(h->index2)[idx] = size + 2;

    /* maintain heap */
    igraph_i_gen2wheap_shift_up(h, size);

    return IGRAPH_SUCCESS;
}

/**
 * Returns the maximum number of elements that the two-way heap can hold. This
 * is also one larger than the maximum allowed index that can be passed to
 * \c igraph_gen2wheap_push_with_index .
 */
igraph_integer_t igraph_gen2wheap_max_size(const igraph_gen2wheap_t *h) {
    return h->max_size;
}

/**
 * Returns a pointer to the largest element in the heap.
 */
const void *igraph_gen2wheap_max(const igraph_gen2wheap_t *h) {
    return ELEM(h, 0);
}

/**
 * Returns the index that was associated to the largest element in the heap
 * when it was pushed to the heap.
 */
igraph_integer_t igraph_gen2wheap_max_index(const igraph_gen2wheap_t *h) {
    return VECTOR(h->index)[0];
}

/**
 * Returns whether the heap contains an element with the given index, even if
 * it was deactivated earlier.
 */
igraph_bool_t igraph_gen2wheap_has_elem(const igraph_gen2wheap_t *h, igraph_integer_t idx) {
    return VECTOR(h->index2)[idx] != 0;
}

/**
 * Returns whether the heap contains an element with the given index \em and it
 * has not been deactivated yet.
 */
igraph_bool_t igraph_gen2wheap_has_active(const igraph_gen2wheap_t *h, igraph_integer_t idx) {
    return VECTOR(h->index2)[idx] > 1;
}

/**
 * Returns a pointer to the item at the given index in the two-way heap.
 */
const void *igraph_gen2wheap_get(const igraph_gen2wheap_t *h, igraph_integer_t idx) {
    igraph_integer_t i = VECTOR(h->index2)[idx] - 2;
    return ELEM(h, i);
}

/**
 * Deletes the largest element from the two-way heap.
 *
 * This function does \em not change the indices associated to the elements
 * that remain in the heap.
 */
void igraph_gen2wheap_delete_max(igraph_gen2wheap_t *h) {
    igraph_integer_t tmpidx = VECTOR(h->index)[0];
    igraph_i_gen2wheap_switch(h, 0, igraph_gen2wheap_size(h) - 1);
    igraph_vector_int_pop_back(&h->index);
    VECTOR(h->index2)[tmpidx] = 0;
    igraph_i_gen2wheap_sink(h, 0);
}

/**
 * Deactivates the largest element from the two-way heap.
 *
 * This function does \em not change the indices associated to the elements
 * that remain in the heap.
 */
void igraph_gen2wheap_deactivate_max(igraph_gen2wheap_t *h) {
    igraph_integer_t tmpidx = VECTOR(h->index)[0];
    igraph_i_gen2wheap_switch(h, 0, igraph_gen2wheap_size(h) - 1);
    igraph_vector_int_pop_back(&h->index);
    VECTOR(h->index2)[tmpidx] = 1;
    igraph_i_gen2wheap_sink(h, 0);
}

/**
 * Modifies the value associated to the given index in the two-way heap.
 */
void igraph_gen2wheap_modify(igraph_gen2wheap_t *h, igraph_integer_t idx, const void *elem) {

    igraph_integer_t pos = VECTOR(h->index2)[idx] - 2;

    memcpy(ELEM(h, pos), elem, h->item_size);
    igraph_i_gen2wheap_sink(h, pos);
    igraph_i_gen2wheap_shift_up(h, pos);
}

/**
 * Checks that the heap is in a consistent state
 */
igraph_error_t igraph_gen2wheap_check(const igraph_gen2wheap_t *h) {
    igraph_integer_t size = igraph_gen2wheap_size(h);
    igraph_integer_t i;
    igraph_bool_t error = false;

    /* Check the heap property */
    for (i = 0; i < size; i++) {
        if (LEFTCHILD(i) >= size) {
            break;
        }
        if (h->cmp(ELEM(h, LEFTCHILD(i)), ELEM(h, i)) > 0) {
            error = true; break;
        }
        if (RIGHTCHILD(i) >= size) {
            break;
        }
        if (h->cmp(ELEM(h, RIGHTCHILD(i)), ELEM(h, i)) > 0) {
            error = true; break;
        }
    }

    if (error) {
        IGRAPH_ERROR("Inconsistent heap.", IGRAPH_EINTERNAL);
    }

    return IGRAPH_SUCCESS;
}
