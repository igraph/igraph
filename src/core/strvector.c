/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_types.h"
#include "igraph_strvector.h"
#include "igraph_memory.h"
#include "igraph_error.h"

#include "internal/hacks.h" /* strdup */
#include "math/safe_intop.h"

#include <string.h>         /* memcpy & co. */
#include <stdlib.h>

/**
 * \section igraph_strvector_t
 * <para>
 * The <type>igraph_strvector_t</type> type is a vector of null-terminated
 * strings. It is used internally for storing graph attribute names as well as
 * string attributes in the C attribute handler.
 * </para>
 *
 * <para>
 * This container automatically manages the memory of its elements.
 * The strings within an <type>igraph_strvector_t</type> should be considered
 * constant, and not modified directly. Functions that add new elements
 * always make copies of the string passed to them.
 * </para>
 *
 * <para>
 * \example examples/simple/igraph_strvector.c
 * </para>
 */

/**
 * \ingroup strvector
 * \function igraph_strvector_init
 * \brief Initializes a string vector.
 *
 * Reserves memory for the string vector, a string vector must be
 * first initialized before calling other functions on it.
 * All elements of the string vector are set to the empty string.
 *
 * \param sv Pointer to an initialized string vector.
 * \param len The (initial) length of the string vector.
 * \return Error code.
 *
 * Time complexity: O(\p len).
 */

igraph_error_t igraph_strvector_init(igraph_strvector_t *sv, igraph_integer_t size) {

    sv->stor_begin = IGRAPH_CALLOC(size, char*);
    IGRAPH_CHECK_OOM(sv->stor_begin, "Cannot initialize string vector.");

    sv->stor_end = sv->stor_begin + size;
    sv->end = sv->stor_end;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_destroy
 * \brief Frees the memory allocated for the string vector.
 *
 * Destroy a string vector. It may be reinitialized with \ref
 * igraph_strvector_init() later.
 * \param sv The string vector.
 *
 * Time complexity: O(l), the total length of the strings, maybe less
 * depending on the memory manager.
 */

void igraph_strvector_destroy(igraph_strvector_t *sv) {
    char **ptr;
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    for (ptr = sv->stor_begin; ptr < sv->end; ptr++) {
        IGRAPH_FREE(*ptr);
    }
    IGRAPH_FREE(sv->stor_begin);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_get
 * \brief Retrieves an element of a string vector.
 *
 * Query an element of a string vector. The returned string must not be modified.
 *
 * \param sv The input string vector.
 * \param idx The index of the element to query.
 *
 * Time complexity: O(1).
 */

const char *igraph_strvector_get(const igraph_strvector_t *sv, igraph_integer_t idx) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    return sv->stor_begin[idx] ? sv->stor_begin[idx] : "";
}

/**
 * \ingroup strvector
 * \function igraph_strvector_set
 * \brief Sets an element of the string vector from a string.
 *
 * The provided \p value is copied into the \p idx position in the
 * string vector.
 *
 * \param sv The string vector.
 * \param idx The position to set.
 * \param value The new value.
 * \return Error code.
 *
 * Time complexity: O(l), the length of the new string. Maybe more,
 * depending on the memory management, if reallocation is needed.
 */

igraph_error_t igraph_strvector_set(igraph_strvector_t *sv, igraph_integer_t idx,
                         const char *value) {
    return igraph_strvector_set_len(sv, idx, value, strlen(value));
}

/**
 * \ingroup strvector
 * \function igraph_strvector_set_len
 * \brief Sets an element of the string vector given a buffer and its size.
 *
 * This is almost the same as \ref igraph_strvector_set, but the new
 * value is not a zero terminated string, but its length is given.
 *
 * \param sv The string vector.
 * \param idx The position to set.
 * \param value The new value.
 * \param len The length of the new value.
 * \return Error code.
 *
 * Time complexity: O(l), the length of the new string. Maybe more,
 * depending on the memory management, if reallocation is needed.
 */
igraph_error_t igraph_strvector_set_len(igraph_strvector_t *sv, igraph_integer_t idx,
                          const char *value, size_t len) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    if (sv->stor_begin[idx] == NULL) {
        sv->stor_begin[idx] = strndup(value, len);
        IGRAPH_CHECK_OOM(sv->stor_begin[idx], "Cannot reserve space for new item in string vector.");
    } else {
        char *tmp = IGRAPH_REALLOC(sv->stor_begin[idx], len + 1, char);
        IGRAPH_CHECK_OOM(tmp, "Cannot reserve space for new item in string vector.");

        sv->stor_begin[idx] = tmp;
        memcpy(sv->stor_begin[idx], value, len * sizeof(char));
        sv->stor_begin[idx][len] = '\0';
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_remove_section
 * \brief Removes a section from a string vector.
 *
 * This function removes the range <code>[from, to)</code> from the string vector.
 *
 * \param sv The string vector.
 * \param from The position of the first element to remove.
 * \param to   The position of the first element \em not to remove.
 */

void igraph_strvector_remove_section(
        igraph_strvector_t *sv, igraph_integer_t from, igraph_integer_t to) {
    igraph_integer_t size = igraph_strvector_size(sv);
    igraph_integer_t i;

    if (from < 0) {
        from = 0;
    }

    if (to > size) {
        to = size;
    }

    if (to > from) {
        for (i = from; i < to; i++) {
            IGRAPH_FREE(sv->stor_begin[i]);
        }

        memmove(sv->stor_begin + from, sv->stor_begin + to,
                sizeof(char*) * (sv->end - sv->stor_begin - to));
        sv->end -= (to - from);
    }
}

/**
 * \ingroup strvector
 * \function igraph_strvector_remove
 * \brief Removes a single element from a string vector.
 *
 * The string will be one shorter.
 * \param sv The string vector.
 * \param elem The index of the element to remove.
 *
 * Time complexity: O(n), the length of the string.
 */

void igraph_strvector_remove(igraph_strvector_t *sv, igraph_integer_t elem) {
    igraph_strvector_remove_section(sv, elem, elem + 1);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_init_copy
 * \brief Initialization by copying.
 *
 * Initializes a string vector by copying another string vector.
 *
 * \param to Pointer to an uninitialized string vector.
 * \param from The other string vector, to be copied.
 * \return Error code.
 *
 * Time complexity: O(l), the total length of the strings in \p from.
 */

igraph_error_t igraph_strvector_init_copy(igraph_strvector_t *to,
                                          const igraph_strvector_t *from) {
    igraph_integer_t from_size = igraph_strvector_size(from);

    to->stor_begin = IGRAPH_CALLOC(from_size, char*);
    IGRAPH_CHECK_OOM(to->stor_begin, "Cannot copy string vector.");

    for (igraph_integer_t i = 0; i < from_size; i++) {
        /* If the string in the 'from' vector is empty, we represent it as NULL.
         * The NULL value was already set by IGRAPH_CALLOC(). */
        if (from->stor_begin[i] == NULL || from->stor_begin[i][0] == '\0') {
            continue;
        }
        to->stor_begin[i] = strdup(from->stor_begin[i]);
        if (to->stor_begin[i] == NULL) {
            /* LCOV_EXCL_START */
            for (igraph_integer_t j = 0; j < i; j++) {
                IGRAPH_FREE(to->stor_begin[j]);
            }
            IGRAPH_FREE(to->stor_begin);
            IGRAPH_ERROR("Cannot copy string vector.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
            /* LCOV_EXCL_STOP */
        }
    }

    to->stor_end = to->stor_begin + from_size;
    to->end = to->stor_end;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_copy
 * \brief Initialization by copying (deprecated alias).
 *
 * \deprecated-by igraph_strvector_init_copy 0.10.0
 */

igraph_error_t igraph_strvector_copy(igraph_strvector_t *to,
                          const igraph_strvector_t *from) {
    return igraph_strvector_init_copy(to, from);
}

/**
 * \function igraph_strvector_append
 * \brief Concatenates two string vectors.
 *
 * Appends the contents of the \p from vector to the \p to vector.
 * If the \p from vector is no longer needed after this operation,
 * use \ref igraph_strvector_merge() for better performance.
 *
 * \param to The first string vector, the result is stored here.
 * \param from The second string vector, it is kept unchanged.
 * \return Error code.
 *
 * \sa \ref igraph_strvector_merge()
 *
 * Time complexity: O(n+l2), n is the number of strings in the new
 * string vector, l2 is the total length of strings in the \p from
 * string vector.
 */

igraph_error_t igraph_strvector_append(igraph_strvector_t *to,
                            const igraph_strvector_t *from) {
    igraph_integer_t len1 = igraph_strvector_size(to), len2 = igraph_strvector_size(from);
    igraph_integer_t newlen;
    igraph_bool_t error = false;
    char *tmp;

    IGRAPH_SAFE_ADD(len1, len2, &newlen);
    IGRAPH_CHECK(igraph_strvector_reserve(to, newlen));

    for (igraph_integer_t i = 0; i < len2; i++) {
        if (from->stor_begin[i] == NULL || from->stor_begin[i][0] == '\0') {
            /* Represent empty strings as NULL. */
            tmp = NULL;
        } else {
            tmp = strdup(from->stor_begin[i]);
            if (tmp == NULL) {
                error = true;
                break;
            }
        }
        *(to->end) = tmp;
        to->end++;
    }

    if (error) {
        igraph_strvector_resize(to, len1); /* always shrinks */
        IGRAPH_ERROR("Cannot append string vector.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_merge
 * \brief Moves the contents of a string vector to the end of another.
 *
 * Transfers the contents of the \p from vector to the end of \p to, clearing
 * \p from in the process. If this operation fails, both vectors are left intact.
 * This function does not copy or reallocate individual strings, therefore it
 * performs better than \ref igraph_strvector_append().
 *
 * \param to   The target vector. The contents of \p from will be appended to it.
 * \param from The source vector. It will be cleared.
 * \return Error code.
 *
 * \sa \ref igraph_strvector_append()
 *
 * Time complexity: O(l2) if \p to has sufficient capacity, O(2*l1+l2) otherwise,
 *   where l1 and l2 are the lengths of \p to and \from respectively.
 */
igraph_error_t igraph_strvector_merge(igraph_strvector_t *to, igraph_strvector_t *from) {
    char **p1, **p2, **pe;
    igraph_integer_t newlen;

    IGRAPH_SAFE_ADD(igraph_strvector_size(to), igraph_strvector_size(from), &newlen);
    IGRAPH_CHECK(igraph_strvector_reserve(to, newlen));

    /* transfer contents of 'from' to 'to */
    for (p1 = to->end, p2 = from->stor_begin, pe = to->stor_begin + newlen;
         p1 < pe; ++p1, ++p2)
    {
        *p1 = *p2;
    }
    to->end = pe;

    /* clear 'from' */
    from->end = from->stor_begin;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_strvector_clear
 * \brief Removes all elements from a string vector.
 *
 * After this operation the string vector will be empty.
 *
 * \param sv The string vector.
 *
 * Time complexity: O(l), the total length of strings, maybe less,
 * depending on the memory manager.
 */

void igraph_strvector_clear(igraph_strvector_t *sv) {
    igraph_integer_t n = igraph_strvector_size(sv);

    for (igraph_integer_t i = 0; i < n; i++) {
        IGRAPH_FREE(sv->stor_begin[i]);
    }
    sv->end = sv->stor_begin;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_resize
 * \brief Resizes a string vector.
 *
 * If the new size is bigger then empty strings are added, if it is
 * smaller then the unneeded elements are removed.
 *
 * \param sv The string vector.
 * \param newsize The new size.
 * \return Error code.
 *
 * Time complexity: O(n), the number of strings if the vector is made
 * bigger, O(l), the total length of the deleted strings if it is made
 * smaller, maybe less, depending on memory management.
 */

igraph_error_t igraph_strvector_resize(igraph_strvector_t *sv, igraph_integer_t newsize) {
    igraph_integer_t toadd = newsize - igraph_strvector_size(sv);
    igraph_integer_t oldsize = igraph_strvector_size(sv);

    if (newsize < oldsize) {
        for (igraph_integer_t i = newsize; i < oldsize; i++) {
            IGRAPH_FREE(sv->stor_begin[i]);
        }
        sv->end = sv->stor_begin + newsize;
    } else if (newsize > oldsize) {
        IGRAPH_CHECK(igraph_strvector_reserve(sv, newsize));
        memset(sv->stor_begin + oldsize, 0, toadd * sizeof(char *));
        sv->end = sv->stor_begin + newsize;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_capacity
 * \brief Returns the capacity of a string vector.
 *
 * \param sv The string vector.
 * \return The capacity of the string vector.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_strvector_capacity(const igraph_strvector_t *sv) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    return sv->stor_end - sv->stor_begin;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_reserve
 * \brief Reserves memory for a string vector.
 *
 * </para><para>
 * \a igraph string vectors are flexible, they can grow and
 * shrink. Growing however occasionally needs the data in the vector to be copied.
 * In order to avoid this, you can call this function to reserve space for
 * future growth of the vector.
 *
 * </para><para>
 * Note that this function does \em not change the size of the
 * string vector. Let us see a small example to clarify things: if you
 * reserve space for 100 strings and the size of your
 * vector was (and still is) 60, then you can surely add additional 40
 * strings to your vector before it will be copied.
 *
 * \param sv The string vector object.
 * \param capacity The new \em allocated size of the string vector.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the new allocated size of the vector.
 */

igraph_error_t igraph_strvector_reserve(igraph_strvector_t *sv, igraph_integer_t capacity) {
    igraph_integer_t current_capacity = igraph_strvector_capacity(sv);
    char **tmp;

    if (capacity <= current_capacity) {
        return IGRAPH_SUCCESS;
    }

    tmp = IGRAPH_REALLOC(sv->stor_begin, capacity, char *);
    IGRAPH_CHECK_OOM(tmp, "Cannot reserve space for new items in string vector.");

    sv->end = tmp + (sv->end - sv->stor_begin);
    sv->stor_begin = tmp;
    sv->stor_end = sv->stor_begin + capacity;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_resize_min
 * \brief Deallocates the unused memory of a string vector.
 *
 * This function attempts to deallocate the unused reserved storage
 * of a string vector. If it succeeds, \ref igraph_strvector_size() and
 * \ref igraph_strvector_capacity() will be the same. The data in the
 * string vector is always preserved, even if deallocation is not successful.
 *
 * \param sv The string vector.
 *
 * Time complexity: Operating system dependent, at most O(n).
 */

void igraph_strvector_resize_min(igraph_strvector_t *sv) {
    igraph_integer_t size;
    char **tmp;
    if (sv->stor_end == sv->end) {
        return;
    }

    size = (sv->end - sv->stor_begin);
    tmp = IGRAPH_REALLOC(sv->stor_begin, size, char *);

    if (tmp != NULL) {
        sv->stor_begin = tmp;
        sv->stor_end = sv->end = sv->stor_begin + size;
    }
}

/**
 * \ingroup strvector
 * \function igraph_strvector_size
 * \brief Returns the size of a string vector.
 *
 * \param sv The string vector.
 * \return The length of the string vector.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_strvector_size(const igraph_strvector_t *sv) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    return sv->end - sv->stor_begin;
}

/**
 * Ensures that the vector has at least one extra slot at the end of its
 * allocated storage area.
 */
static igraph_error_t igraph_i_strvector_expand_if_full(igraph_strvector_t *sv) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    if (sv->stor_end == sv->end) {
        igraph_integer_t old_size = igraph_strvector_size(sv);
        igraph_integer_t new_size = old_size < IGRAPH_INTEGER_MAX/2 ? old_size * 2 : IGRAPH_INTEGER_MAX;
        if (old_size == IGRAPH_INTEGER_MAX) {
            IGRAPH_ERROR("Cannot add new item to string vector, already at maximum size.", IGRAPH_EOVERFLOW);
        }
        if (new_size == 0) {
            new_size = 1;
        }
        IGRAPH_CHECK(igraph_strvector_reserve(sv, new_size));
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_push_back
 * \brief Adds an element to the back of a string vector.
 *
 * \param sv The string vector.
 * \param value The string to add; it will be copied.
 * \return Error code.
 *
 * Time complexity: O(n+l), n is the total number of strings, l is the
 * length of the new string.
 */

igraph_error_t igraph_strvector_push_back(igraph_strvector_t *sv, const char *value) {
    IGRAPH_CHECK(igraph_i_strvector_expand_if_full(sv));
    char *tmp = strdup(value);
    IGRAPH_CHECK_OOM(tmp, "Cannot push new string to string vector.");
    *sv->end = tmp;
    sv->end++;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_push_back_len
 * \brief Adds a string of the given length to the back of a string vector.
 *
 * \param sv The string vector.
 * \param value The start of the string to add. At most \p len characters will be copied.
 * \param len The length of the string.
 * \return Error code.
 *
 * Time complexity: O(n+l), n is the total number of strings, l is the
 * length of the new string.
 */

igraph_error_t igraph_strvector_push_back_len(
        igraph_strvector_t *sv,
        const char *value, igraph_integer_t len) {

    IGRAPH_CHECK(igraph_i_strvector_expand_if_full(sv));
    char *tmp = strndup(value, len);
    if (! tmp) {
        IGRAPH_ERROR("Cannot add string to string vector.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    *sv->end = tmp;
    sv->end++;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_add
 * \brief Adds an element to the back of a string vector (deprecated alias).
 *
 * \deprecated-by igraph_strvector_push_back 0.10.0
 */

igraph_error_t igraph_strvector_add(igraph_strvector_t *sv, const char *value) {
    return igraph_strvector_push_back(sv, value);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_set2
 * \brief Sets an element of the string vector given a buffer and its size (deprecated alias).
 *
 * \deprecated-by igraph_strvector_set_len 0.10.0
 */

igraph_error_t igraph_strvector_set2(
    igraph_strvector_t *sv, igraph_integer_t idx, const char *value, size_t len
) {
    return igraph_strvector_set_len(sv, idx, value, len);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_print
 * \brief Prints a string vector.
 *
 * \param sv The string vector.
 * \param file The file to write to.
 * \param sep The separator to print between strings.
 * \return Error code.
 */
igraph_error_t igraph_strvector_print(const igraph_strvector_t *sv, FILE *file,
                           const char *sep) {

    igraph_integer_t n = igraph_strvector_size(sv);
    if (n != 0) {
        fprintf(file, "%s", igraph_strvector_get(sv, 0));
    }
    for (igraph_integer_t i = 1; i < n; i++) {
        fprintf(file, "%s%s", sep, igraph_strvector_get(sv, i));
    }
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_index
 * \brief Takes elements at given positions from a string vector.
 *
 * \param sv   The string vector.
 * \param newv An initialized string vector, it will be resized as needed.
 * \param idx  An integer vector of indices to take from \p sv.
 * \return Error code.
 */
igraph_error_t igraph_strvector_index(const igraph_strvector_t *sv,
                           igraph_strvector_t *newv,
                           const igraph_vector_int_t *idx) {

    igraph_integer_t newlen = igraph_vector_int_size(idx);
    IGRAPH_CHECK(igraph_strvector_resize(newv, newlen));

    for (igraph_integer_t i = 0; i < newlen; i++) {
        igraph_integer_t j = VECTOR(*idx)[i];
        const char *str = igraph_strvector_get(sv, j);
        IGRAPH_CHECK(igraph_strvector_set(newv, i, str));
    }

    return IGRAPH_SUCCESS;
}
