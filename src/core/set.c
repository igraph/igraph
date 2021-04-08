/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"

#include "core/set.h"

#include <string.h>     /* memmove */

#define SET(s) ((s).stor_begin)

/**
 * \ingroup set
 * \function igraph_set_init
 * \brief Initializes a set.
 *
 * \param set pointer to the set to be initialized
 * \param size the expected number of elements in the set
 *
 * \return error code:
 *       \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the expected size of the set.
 */
int igraph_set_init(igraph_set_t *set, int long size) {
    long int alloc_size;
    
    if (size < 0) {
        size = 0;
    }
    alloc_size = size > 0 ? size : 1;
    set->stor_begin = IGRAPH_CALLOC(alloc_size, igraph_integer_t);
    set->stor_end = set->stor_begin + alloc_size;
    set->end = set->stor_begin;

    return 0;
}

/**
 * \ingroup set
 * \function igraph_set_destroy
 * \brief Destroys a set object.
 *
 * \param set pointer to the set to be destroyed
 *
 * Time complexity: operating system dependent.
 */
void igraph_set_destroy(igraph_set_t* set) {
    IGRAPH_ASSERT(set != 0);
    if (set->stor_begin != 0) {
        IGRAPH_FREE(set->stor_begin);
        set->stor_begin = NULL;
    }
}

/**
 * \ingroup set
 * \function igraph_set_inited
 * \brief Determines whether a set is initialized or not.
 *
 * This function checks whether the internal storage for the members of the
 * set has been allocated or not, and it assumes that the pointer for the
 * internal storage area contains \c NULL if the area is not initialized yet.
 * This only applies if you have allocated an array of sets with \c IGRAPH_CALLOC or
 * if you used the \c IGRAPH_SET_NULL constant to initialize the set.
 *
 * \param set The set object.
 *
 * Time complexity: O(1)
 */
igraph_bool_t igraph_set_inited(igraph_set_t* set) {
    return (set->stor_begin != 0);
}

/**
 * \ingroup set
 * \function igraph_set_reserve
 * \brief Reserve memory for a set.
 *
 * \param set The set object.
 * \param size the new \em allocated size of the set.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the new allocated size of the set.
 */
int igraph_set_reserve(igraph_set_t* set, long int size) {
    long int actual_size = igraph_set_size(set);
    igraph_integer_t *tmp;
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    if (size <= actual_size) {
        return 0;
    }

    tmp = IGRAPH_REALLOC(set->stor_begin, (size_t) size, igraph_integer_t);
    if (tmp == 0) {
        IGRAPH_ERROR("cannot reserve space for set", IGRAPH_ENOMEM);
    }
    set->stor_begin = tmp;
    set->stor_end = set->stor_begin + size;
    set->end = set->stor_begin + actual_size;

    return 0;
}

/**
 * \ingroup set
 * \function igraph_set_empty
 * \brief Decides whether the size of the set is zero.
 *
 * \param set The set object.
 * \return Non-zero number if the size of the set is not zero and
 *         zero otherwise.
 *
 * Time complexity: O(1).
 */
igraph_bool_t igraph_set_empty(const igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    return set->stor_begin == set->end;
}

/**
 * \ingroup set
 * \function igraph_set_clear
 * \brief Removes all elements from a set.
 *
 * </para><para>
 * This function simply sets the size of the set to zero, it does
 * not free any allocated memory. For that you have to call
 * \ref igraph_set_destroy().
 * \param v The set object.
 *
 * Time complexity: O(1).
 */
void igraph_set_clear(igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    set->end = set->stor_begin;
}


/**
 * \ingroup set
 * \function igraph_set_size
 * \brief Gives the size (=length) of the set.
 *
 * \param v The set object
 * \return The size of the set.
 *
 * Time complexity: O(1).
 */

long int igraph_set_size(const igraph_set_t* set) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    return set->end - set->stor_begin;
}


/**
 * \ingroup set
 * \function igraph_set_add
 * \brief Adds an element to the set.
 *
 * \param set The set object.
 * \param e The element to be added.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
int igraph_set_add(igraph_set_t* set, igraph_integer_t e) {
    long int left, right, middle;
    long int size;
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);

    size = igraph_set_size(set);

    /* search where to insert the new element */
    left = 0;
    right = size - 1;
    while (left < right - 1) {
        middle = (left + right) / 2;
        if (SET(*set)[middle] > e) {
            right = middle;
        } else if (SET(*set)[middle] < e) {
            left = middle;
        } else {
            left = middle;
            break;
        }
    }

    if (right >= 0 && SET(*set)[left] != e && SET(*set)[right] == e) {
        left = right;
    }

    while (left < size && set->stor_begin[left] < e) {
        left++;
    }
    if (left >= size || set->stor_begin[left] != e) {
        /* full, allocate more storage */
        if (set->stor_end == set->end) {
            long int new_size = size * 2;
            if (new_size == 0) {
                new_size = 1;
            }
            IGRAPH_CHECK(igraph_set_reserve(set, new_size));
        }

        /* Element should be inserted at position 'left' */
        if (left < size)
            memmove(set->stor_begin + left + 1, set->stor_begin + left,
                    (size_t) (size - left)*sizeof(set->stor_begin[0]));

        set->stor_begin[left] = e;
        set->end += 1;
    }

    return 0;
}

/**
 * \ingroup set
 * \function igraph_set_contains
 * \brief Checks whether a given element is in the set or not.
 *
 * \param set The set object.
 * \param e The element being sought.
 * \return Positive integer (true) if \p e is found, zero (false) otherwise.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
igraph_bool_t igraph_set_contains(igraph_set_t* set, igraph_integer_t e) {
    long int left, right, middle;

    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);

    left = 0;
    right = igraph_set_size(set) - 1;

    if (right == -1) {
        return 0;    /* the set is empty */
    }

    /* search for the new element */
    while (left < right - 1) {
        middle = (left + right) / 2;
        if (SET(*set)[middle] > e) {
            right = middle;
        } else if (SET(*set)[middle] < e) {
            left = middle;
        } else {
            return 1;
        }
    }

    return SET(*set)[left] == e || SET(*set)[right] == e;
}

/**
 * \ingroup set
 * \function igraph_set_iterate
 * \brief Iterates through the element to the set.
 *
 * Elements are returned in an arbitrary order.
 *
 * \param set The set object.
 * \param state Internal state of the iteration.
 *   This should be a pointer to a \c long variable
 *   which must be zero for the first invocation.
 *   The object should not be adjusted and its value should
 *   not be used for anything during the iteration.
 * \param element The next element or \c NULL (if the iteration
 *   has ended) is returned here.
 *
 * \return Nonzero if there are more elements, zero otherwise.
 */
igraph_bool_t igraph_set_iterate(igraph_set_t* set, long int* state,
                                 igraph_integer_t* element) {
    IGRAPH_ASSERT(set != 0);
    IGRAPH_ASSERT(set->stor_begin != 0);
    IGRAPH_ASSERT(state != 0);
    IGRAPH_ASSERT(element != 0);

    if (*state < igraph_set_size(set)) {
        *element = set->stor_begin[*state];
        *state = *state + 1;
        return 1;
    } else {
        *element = 0;
        return 0;
    }
}

