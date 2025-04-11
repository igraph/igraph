/*
   IGraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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

#include "core/set.h"

#include <string.h>     /* memmove */

#define SET(s) ((s).stor_begin)

/**
 * \ingroup set
 * \function igraph_set_init
 * \brief Initializes a set.
 *
 * Initializes an empty set (with zero elements). Allocates memory for
 * the requested capacity. No re-allocation will be necessary until the
 * number of elements exceeds this initial capacity.
 *
 * \param set Pointer to the set to be initialized.
 * \param capacity The expected number of elements in the set.
 *
 * \return error code:
 *       \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the expected size of the set.
 */
igraph_error_t igraph_set_init(igraph_set_t *set, igraph_integer_t capacity) {
    igraph_integer_t alloc_size;

    IGRAPH_ASSERT(capacity >= 0);
    alloc_size = capacity > 0 ? capacity : 1;
    set->stor_begin = IGRAPH_CALLOC(alloc_size, igraph_integer_t);
    if (! set->stor_begin) {
        IGRAPH_ERROR("Cannot initialize set.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    set->stor_end = set->stor_begin + alloc_size;
    set->end = set->stor_begin;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup set
 * \function igraph_set_destroy
 * \brief Destroys a set object.
 *
 * \param set Pointer to the set to be destroyed.
 *
 * Time complexity: operating system dependent.
 */
void igraph_set_destroy(igraph_set_t *set) {
    IGRAPH_ASSERT(set != NULL);
    if (set->stor_begin != NULL) {
        IGRAPH_FREE(set->stor_begin); /* sets to NULL */
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
igraph_bool_t igraph_set_inited(igraph_set_t *set) {
    return (set->stor_begin != NULL);
}

/**
 * \ingroup set
 * \function igraph_set_reserve
 * \brief Reserves memory for a set.
 *
 * \param set The set object.
 * \param capacity the new \em allocated capacity of the set.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the new allocated size of the set.
 */
igraph_error_t igraph_set_reserve(igraph_set_t *set, igraph_integer_t capacity) {
    igraph_integer_t actual_size = igraph_set_size(set);
    igraph_integer_t *tmp;
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    if (capacity <= actual_size) {
        return IGRAPH_SUCCESS;
    }

    tmp = IGRAPH_REALLOC(set->stor_begin, capacity, igraph_integer_t);
    IGRAPH_CHECK_OOM(tmp, "Cannot reserve space for set.");

    set->stor_begin = tmp;
    set->stor_end = set->stor_begin + capacity;
    set->end = set->stor_begin + actual_size;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup set
 * \function igraph_set_empty
 * \brief Decides whether the size of the set is zero.
 *
 * \param set The set object.
 * \return True if the size of the set is not zero and false otherwise.
 *
 * Time complexity: O(1).
 */
igraph_bool_t igraph_set_empty(const igraph_set_t *set) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    return set->stor_begin == set->end;
}

/**
 * \ingroup set
 * \function igraph_set_clear
 * \brief Removes all elements from the set.
 *
 * This function simply sets the size of the set to zero, it does
 * not free any allocated memory. For that you have to call
 *
 * \ref igraph_set_destroy().
 *
 * \param set The set object.
 *
 * Time complexity: O(1).
 */
void igraph_set_clear(igraph_set_t *set) {
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);
    set->end = set->stor_begin;
}


/**
 * \ingroup set
 * \function igraph_set_size
 * \brief Gives the size of the set.
 *
 * The number of elements in the set.
 *
 * \param set The set object
 * \return The size of the set.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_set_size(const igraph_set_t *set) {
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
igraph_error_t igraph_set_add(igraph_set_t *set, igraph_integer_t e) {
    igraph_integer_t left, right, middle;
    igraph_integer_t size;
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
            igraph_integer_t new_size = size < IGRAPH_INTEGER_MAX/2 ? size * 2 : IGRAPH_INTEGER_MAX;
            if (size == IGRAPH_INTEGER_MAX) {
                IGRAPH_ERROR("Cannot add to set, already at maximum size.", IGRAPH_EOVERFLOW);
            }
            if (new_size == 0) {
                new_size = 1;
            }
            IGRAPH_CHECK(igraph_set_reserve(set, new_size));
        }

        /* Element should be inserted at position 'left' */
        if (left < size)
            memmove(set->stor_begin + left + 1, set->stor_begin + left,
                    (size - left) * sizeof(set->stor_begin[0]));

        set->stor_begin[left] = e;
        set->end += 1;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup set
 * \function igraph_set_remove
 * \brief Adds an element to the set.
 *
 * If element not present, does nothing. Does not change capacity!
 *
 * \param set The set object.
 * \param e The element to be removed.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
void igraph_set_remove(igraph_set_t *set, igraph_integer_t e) {
    igraph_integer_t left, right, middle;
    igraph_integer_t size;
    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);

    size = igraph_set_size(set);
    if (size == 0) return;

    /* find the element to remove */
    // TODO: extract binsearch to own function for add and remove functions
    // Also why is there an in-house binsearch implementation...
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

    // make sure left points at e
    if (right >= 0 && SET(*set)[right] == e) {
        left = right;
    }

    // if it still doesn't then element not found, there's nothing to do, return
    if (SET(*set)[left] != e) return;

    // if the element we found is not the last one, shift all following elements
    // to the left by one.
    if (left + 1 < size) {
        memmove(
            set->stor_begin + left,
            set->stor_begin + left + 1,
            (size - left - 1) * sizeof(set->stor_begin[0])
        );
    }

    // delete the last element and return
    set->end -= 1;
}

/**
 * \ingroup set
 * \function igraph_set_remove
 * \brief Remove a set of elements from a set.
 *
 * Modifies \c set in-place. Computes \c set-to_remove
 *
 * \param set The set object.
 * \param to_remove The element to be removed.
 *
 * Time complexity: O(n+m), n and m are the number of elements in \p set and
 * \p to_remove respectively.
 */
void igraph_set_difference(igraph_set_t *set, igraph_set_t *to_remove) {
    if (igraph_set_empty(set) || igraph_set_empty(to_remove)) return;  // nothing to do

    igraph_integer_t a_read = 0, a_write = 0;  // read and write indices into set
    igraph_integer_t b_read = 0;  // read index into to_remove
    igraph_integer_t a_size = igraph_set_size(set);
    igraph_integer_t b_size = igraph_set_size(to_remove);
    igraph_integer_t unread, removed = 0;

    igraph_integer_t *A = set->stor_begin;
    igraph_integer_t *B = to_remove->stor_begin;

    while (true) {
        // catch b_read up to current a_read state
        while (b_read < b_size && B[b_read] < A[a_read]) {
            b_read++;
        }

        // if ran out of B elements, done
        if (b_read == b_size) {
            // there may still be unread A entries, they are all larger than all
            // element of B so they should all be kept; memmove from read to write.
            unread = a_size - a_read;
            if (unread > 0) {
                memmove(
                    set->stor_begin+a_write,
                    set->stor_begin+a_read,
                    unread * sizeof(set->stor_begin[0])
                );
            }
            break;
        }

        // catch a_read up to current b_read state
        while (a_read < a_size && A[a_read] < B[b_read]) {
            // the elements seen in A are not in B so they are safe to write
            A[a_write] = A[a_read];
            a_read++;
            a_write++; // can't overshoot, since a_write <= a_read always holds
        }

        // if ran out of A elements, done
        if (a_read == a_size) {
            break;
        }

        // check if A element in B; if it is, it should be removed
        // that means move read index ahead, but not write index
        if (A[a_read] == B[b_read]) {
            removed++;
            a_read++;
            b_read++;
        }
    }

    // Now a_write points to just after the last element in the new set A\B
    set->end -= removed;
}

/**
 * \ingroup set
 * \function igraph_set_contains
 * \brief Checks whether a given element is in the set or not.
 *
 * \param set The set object.
 * \param e The element being sought.
 * \return True if \p e is found, false otherwise.
 *
 * Time complexity: O(log(n)), n is the number of elements in \p set.
 */
igraph_bool_t igraph_set_contains(const igraph_set_t *set, igraph_integer_t e) {
    igraph_integer_t left, right, middle;

    IGRAPH_ASSERT(set != NULL);
    IGRAPH_ASSERT(set->stor_begin != NULL);

    left = 0;
    right = igraph_set_size(set) - 1;

    if (right == -1) {
        return false;    /* the set is empty */
    }

    /* search for the new element */
    while (left < right - 1) {
        middle = (left + right) / 2;
        if (SET(*set)[middle] > e) {
            right = middle;
        } else if (SET(*set)[middle] < e) {
            left = middle;
        } else {
            return true;
        }
    }

    return SET(*set)[left] == e || SET(*set)[right] == e;
}

/**
 * \ingroup set
 * \function igraph_set_iterate
 * \brief Iterates through the element of the set.
 *
 * Elements are returned in an arbitrary order.
 *
 * \param set The set object.
 * \param state Internal state of the iteration.
 *   This should be a pointer to a \type igraph_integer_t variable
 *   which must be zero for the first invocation.
 *   The object must not be adjusted and its value should
 *   not be used for anything during the iteration.
 * \param element The next element or 0 (if the iteration
 *   has ended) is returned here.
 *
 * \return True if there are more elements, false otherwise.
 */
igraph_bool_t igraph_set_iterate(const igraph_set_t *set, igraph_integer_t *state,
                                 igraph_integer_t *element) {
    IGRAPH_ASSERT(set != 0);
    IGRAPH_ASSERT(set->stor_begin != 0);
    IGRAPH_ASSERT(state != 0);
    IGRAPH_ASSERT(element != 0);

    if (*state < igraph_set_size(set)) {
        *element = set->stor_begin[*state];
        *state = *state + 1;
        return true;
    } else {
        *element = 0;
        return false;
    }
}
