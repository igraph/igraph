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

#include <string.h>         /* memcpy & co. */
#include <stdlib.h>

#include "igraph_types.h"
#include "igraph_strvector.h"
#include "igraph_memory.h"
#include "igraph_error.h"

#include "internal/hacks.h" /* strdup */

/**
 * \section igraph_strvector_t
 * <para>The <type>igraph_strvector_t</type> type is a vector of strings.
 * The current implementation is very simple and not too efficient. It
 * works fine for not too many strings, e.g. the list of attribute
 * names is returned in a string vector by \ref
 * igraph_cattribute_list(). Do not expect great performance from this
 * type.</para>
 *
 * <para>
 * \example examples/simple/igraph_strvector.c
 * </para>
 */

/**
 * \ingroup strvector
 * \function igraph_strvector_init
 * \brief Initialize
 *
 * Reserves memory for the string vector, a string vector must be
 * first initialized before calling other functions on it.
 * All elements of the string vector are set to the empty string.
 * \param sv Pointer to an initialized string vector.
 * \param len The (initial) length of the string vector.
 * \return Error code.
 *
 * Time complexity: O(\p len).
 */

igraph_error_t igraph_strvector_init(igraph_strvector_t *sv, igraph_integer_t size) {
    igraph_integer_t i;
    sv->stor_begin = IGRAPH_CALLOC(size, char*);
    if (sv->stor_begin == NULL) {
        IGRAPH_ERROR("Strvector init failed.", IGRAPH_ENOMEM);
    }
    for (i = 0; i < size; i++) {
        sv->stor_begin[i] = IGRAPH_CALLOC(1, char);
        if (sv->stor_begin[i] == NULL) {
            igraph_strvector_destroy(sv);
            IGRAPH_ERROR("Strvector init failed.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
        sv->stor_begin[i][0] = '\0';
    }
    sv->stor_end = sv->stor_begin + size;
    sv->end = sv->stor_end;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_destroy
 * \brief Free allocated memory
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
 * \brief Indexing
 *
 * Query an element of a string vector. See also the \ref STR macro
 * for an easier way.
 * \param sv The input string vector.
 * \param idx The index of the element to query.
 *
 * Time complexity: O(1).
 */

char* igraph_strvector_get(const igraph_strvector_t *sv, igraph_integer_t idx) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    IGRAPH_ASSERT(sv->stor_begin[idx] != NULL);
    return sv->stor_begin[idx];
}

/**
 * \ingroup strvector
 * \function igraph_strvector_set
 * \brief Set an element
 *
 * The provided \p value is copied into the \p idx position in the
 * string vector.
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
    size_t value_len;
    char *tmp;

    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    IGRAPH_ASSERT(sv->stor_begin[idx] != NULL);

    value_len = strlen(value);
    tmp = IGRAPH_REALLOC(sv->stor_begin[idx], value_len + 1, char);
    if (tmp == NULL) {
        IGRAPH_ERROR("Strvector set failed.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    sv->stor_begin[idx] = tmp;
    strcpy(sv->stor_begin[idx], value);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_set2
 * \brief Sets an element.
 *
 * This is almost the same as \ref igraph_strvector_set, but the new
 * value is not a zero terminated string, but its length is given.
 * \param sv The string vector.
 * \param idx The position to set.
 * \param value The new value.
 * \param len The length of the new value.
 * \return Error code.
 *
 * Time complexity: O(l), the length of the new string. Maybe more,
 * depending on the memory management, if reallocation is needed.
 */
igraph_error_t igraph_strvector_set2(igraph_strvector_t *sv, igraph_integer_t idx,
                          const char *value, size_t len) {
    char *tmp;

    if (idx < 0 || idx >= sv->stor_end - sv->stor_begin) {
        IGRAPH_ERROR("String vector index out of bounds.", IGRAPH_EINVAL);
    }
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    IGRAPH_ASSERT(sv->stor_begin[idx] != NULL);

    tmp = IGRAPH_REALLOC(sv->stor_begin[idx], len + 1, char);
    if (tmp == NULL) {
        IGRAPH_ERROR("Strvector set failed.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    sv->stor_begin[idx] = tmp;
    memcpy(sv->stor_begin[idx], value, len * sizeof(char));
    sv->stor_begin[idx][len] = '\0';

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_remove_section
 * \brief Removes a section from a string vector.
 */

void igraph_strvector_remove_section(
        igraph_strvector_t *sv, igraph_integer_t from, igraph_integer_t to) {
    igraph_integer_t i;

    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    for (i = from; i < to; i++) {
        IGRAPH_FREE(sv->stor_begin[i]);
    }
    for (char **p = sv->stor_begin; p < sv->end - to; p++) {
        p[from] = p[to];
    }

    sv->end -= (to - from);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_remove
 * \brief Removes a single element from a string vector.
 *
 * The string will be one shorter.
 * \param v The string vector.
 * \param elem The index of the element to remove.
 *
 * Time complexity: O(n), the length of the string.
 */

void igraph_strvector_remove(igraph_strvector_t *sv, igraph_integer_t elem) {
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    igraph_strvector_remove_section(sv, elem, elem + 1);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_copy
 * \brief Initialization by copying.
 *
 * Initializes a string vector by copying another string vector.
 * \param to Pointer to an uninitialized string vector.
 * \param from The other string vector, to be copied.
 * \return Error code.
 *
 * Time complexity: O(l), the total length of the strings in \p from.
 */

igraph_error_t igraph_strvector_copy(igraph_strvector_t *to,
                          const igraph_strvector_t *from) {
    igraph_integer_t i;
    IGRAPH_ASSERT(from != NULL);
    /*   IGRAPH_ASSERT(from->stor_begin != 0); */
    to->stor_begin = IGRAPH_CALLOC(igraph_strvector_size(from), char*);
    if (to->stor_begin == NULL) {
        IGRAPH_ERROR("Cannot copy string vector.", IGRAPH_ENOMEM);
    }
    to->stor_end = to->stor_begin + igraph_strvector_size(from);
    to->end = to->stor_end;

    for (i = 0; i < igraph_strvector_size(from); i++) {
        to->stor_begin[i] = strdup(igraph_strvector_get(from, i));
        if (to->stor_begin[i] == NULL) {
            igraph_strvector_destroy(to);
            IGRAPH_ERROR("Cannot copy string vector.", IGRAPH_ENOMEM);
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_strvector_append
 * Concatenate two string vectors.
 *
 * \param to The first string vector, the result is stored here.
 * \param from The second string vector, it is kept unchanged.
 * \return Error code.
 *
 * Time complexity: O(n+l2), n is the number of strings in the new
 * string vector, l2 is the total length of strings in the \p from
 * string vector.
 */

igraph_error_t igraph_strvector_append(igraph_strvector_t *to,
                            const igraph_strvector_t *from) {
    igraph_integer_t len1 = igraph_strvector_size(to), len2 = igraph_strvector_size(from);
    igraph_integer_t i;
    igraph_bool_t error = 0;
    IGRAPH_CHECK(igraph_strvector_resize(to, len1 + len2));
    for (i = 0; i < len2; i++) {
        if (from->stor_begin[i][0] != '\0') {
            IGRAPH_FREE(to->stor_begin[len1 + i]);
            to->stor_begin[len1 + i] = strdup(from->stor_begin[i]);
            if (!to->stor_begin[len1 + i]) {
                error = 1;
                break;
            }
        }
    }
    if (error) {
        igraph_strvector_resize(to, len1);
        IGRAPH_ERROR("Cannot append string vector.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_strvector_clear
 * Remove all elements
 *
 * After this operation the string vector will be empty.
 * \param sv The string vector.
 *
 * Time complexity: O(l), the total length of strings, maybe less,
 * depending on the memory manager.
 */

void igraph_strvector_clear(igraph_strvector_t *sv) {
    igraph_integer_t i, n = igraph_strvector_size(sv);

    for (i = 0; i < n; i++) {
        IGRAPH_FREE(sv->stor_begin[i]);
    }
    sv->end = sv->stor_begin;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_resize
 * \brief Resize
 *
 * If the new size is bigger then empty strings are added, if it is
 * smaller then the unneeded elements are removed.
 * \param v The string vector.
 * \param newsize The new size.
 * \return Error code.
 *
 * Time complexity: O(n), the number of strings if the vector is made
 * bigger, O(l), the total length of the deleted strings if it is made
 * smaller, maybe less, depending on memory management.
 */

igraph_error_t igraph_strvector_resize(igraph_strvector_t *sv, igraph_integer_t newsize) {
    igraph_integer_t toadd = newsize - igraph_strvector_size(sv), i, j;
    igraph_integer_t oldsize = igraph_strvector_size(sv);

    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);
    if (newsize < oldsize) {
        for (i = newsize; i < oldsize; i++) {
            IGRAPH_FREE(sv->stor_begin[i]);
        }
        sv->end = sv->stor_begin + newsize;
    } else if (newsize > oldsize) {
        IGRAPH_CHECK(igraph_strvector_reserve(sv, newsize));

        for (i = 0; i < toadd; i++) {
            sv->stor_begin[oldsize + i] = IGRAPH_CALLOC(1, char);
            if (sv->stor_begin[oldsize + i] == 0) {
                for (j = 0; j < i; j++) {
                    if (sv->stor_begin[oldsize + i] != 0) {
                        IGRAPH_FREE(sv->stor_begin[oldsize + i]);
                    }
                }
                IGRAPH_ERROR("Cannot resize string vector.", IGRAPH_ENOMEM);
            }
            sv->stor_begin[oldsize + i][0] = '\0';
        }
        sv->end = sv->stor_begin + newsize;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_capacity
 * \brief Gives the capacity of a string vector.
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
 * \function igraph_strvector_size
 * \brief Gives the size of a string vector.
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

static igraph_error_t add_one_to_size(igraph_strvector_t *sv, igraph_integer_t old_size) {
    if (sv->end == sv->stor_end) {
        igraph_integer_t new_size;
        new_size = old_size < IGRAPH_INTEGER_MAX/2 ? old_size * 2 : IGRAPH_INTEGER_MAX;
        if (old_size == IGRAPH_INTEGER_MAX) {
            IGRAPH_ERROR("Cannot add to strvector, already at maximum size.", IGRAPH_EOVERFLOW);
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
 * \param v The string vector.
 * \param value The string to add; it will be copied.
 * \return Error code.
 *
 * Time complexity: O(n+l), n is the total number of strings, l is the
 * length of the new string.
 */

igraph_error_t igraph_strvector_push_back(igraph_strvector_t *sv, const char *value) {
    igraph_integer_t old_size;

    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    old_size = igraph_strvector_size(sv);
    IGRAPH_CHECK(add_one_to_size(sv, old_size));
    sv->stor_begin[old_size] = strdup(value);
    if (sv->stor_begin[old_size] == NULL) {
        IGRAPH_ERROR("Cannot add string to string vector.", IGRAPH_ENOMEM);
    }
    sv->end = sv->stor_begin + old_size + 1;

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_push_back_len
 * \brief Adds string of given length to the back of a string vector.
 *
 * \param v The string vector.
 * \param value The start of the string to add. At most \p len characters will be copied.
 * \param len The length of the string.
 * \return Error code.
 *
 * Time complexity: O(n+l), n is the total number of strings, l is the
 * length of the new string.
 */

igraph_error_t igraph_strvector_push_back_len(igraph_strvector_t *sv, const char *value,
        igraph_integer_t len) {
    igraph_integer_t old_size;

    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    old_size = igraph_strvector_size(sv);
    IGRAPH_CHECK(add_one_to_size(sv, old_size));
    sv->stor_begin[old_size] = strndup(value, len);
    if (sv->stor_begin[old_size] == NULL) {
        IGRAPH_ERROR("Cannot add string to string vector.", IGRAPH_ENOMEM);
    }
    sv->end = sv->stor_begin + old_size + 1;

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
 * \function igraph_strvector_permdelete
 * \brief Removes elements from a string vector (for internal use)
 */

void igraph_strvector_permdelete(igraph_strvector_t *sv, const igraph_vector_int_t *index,
                                 igraph_integer_t nremove) {
    igraph_integer_t i;
    char **tmp;
    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    for (i = 0; i < igraph_strvector_size(sv); i++) {
        if (VECTOR(*index)[i] != 0) {
            sv->stor_begin[VECTOR(*index)[i] - 1 ] = sv->stor_begin[i];
        } else {
            IGRAPH_FREE(sv->stor_begin[i]);
        }
    }
    /* Try to make it shorter */
    tmp = IGRAPH_REALLOC(sv->stor_begin, igraph_strvector_size(sv) - nremove ? (igraph_strvector_size(sv) - nremove) : 1, char*);
    if (tmp != NULL) {
        sv->stor_begin = tmp;
    }
    sv->stor_end -= nremove;
    sv->end -= nremove;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_print
 * \brief Prints a string vector.
 *
 * \param v The string vector.
 * \param file The file to write to.
 * \param sep The separator to print between strings.
 * \return Error code.
 */
igraph_error_t igraph_strvector_print(const igraph_strvector_t *sv, FILE *file,
                           const char *sep) {

    igraph_integer_t i, n = igraph_strvector_size(sv);
    if (n != 0) {
        fprintf(file, "%s", STR(*sv, 0));
    }
    for (i = 1; i < n; i++) {
        fprintf(file, "%s%s", sep, STR(*sv, i));
    }
    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_strvector_index(const igraph_strvector_t *sv,
                           igraph_strvector_t *newv,
                           const igraph_vector_int_t *idx) {

    igraph_integer_t i, newlen = igraph_vector_int_size(idx);
    IGRAPH_CHECK(igraph_strvector_resize(newv, newlen));

    for (i = 0; i < newlen; i++) {
        igraph_integer_t j = VECTOR(*idx)[i];
        char *str = igraph_strvector_get(sv, j);
        igraph_strvector_set(newv, i, str);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_reserve
 * \brief Reserves memory for a string vector.
 *
 * </para><para>
 * \a igraph string vectors are flexible, they can grow and
 * shrink. Growing
 * however occasionally needs the data in the vector to be copied.
 * In order to avoid this, you can call this function to reserve space for
 * future growth of the vector.
 *
 * </para><para>
 * Note that this function does \em not change the size of the
 * string vector. Let us see a small example to clarify things: if you
 * reserve space for 100 strings and the size of your
 * vector was (and still is) 60, then you can surely add additional 40
 * strings to your vector before it will be copied.
 * \param sv The string vector object.
 * \param capacity The new \em allocated size of the string vector.
 * \return Error code:
 *         \c IGRAPH_ENOMEM if there is not enough memory.
 *
 * Time complexity: operating system dependent, should be around
 * O(n), n is the new allocated size of the vector.
 */

igraph_error_t igraph_strvector_reserve(igraph_strvector_t *sv, igraph_integer_t capacity) {
    igraph_integer_t current_capacity;
    char **tmp;

    IGRAPH_ASSERT(sv != NULL);
    IGRAPH_ASSERT(sv->stor_begin != NULL);

    current_capacity = igraph_strvector_capacity(sv);

    if (capacity <= current_capacity) {
        return IGRAPH_SUCCESS;
    }

    tmp = IGRAPH_REALLOC(sv->stor_begin, capacity, char *);
    if (tmp == 0) {
        IGRAPH_ERROR("Cannot reserve space for vector.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }

    sv->end = tmp + (sv->end - sv->stor_begin);
    sv->stor_begin = tmp;
    sv->stor_end = sv->stor_begin + capacity;

    return IGRAPH_SUCCESS;
}
