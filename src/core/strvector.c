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

#include <string.h>         /* memcpy & co. */
#include <stdlib.h>

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

int igraph_strvector_init(igraph_strvector_t *sv, long int len) {
    long int i;
    sv->data = IGRAPH_CALLOC(len, char*);
    if (sv->data == 0) {
        IGRAPH_ERROR("strvector init failed", IGRAPH_ENOMEM);
    }
    for (i = 0; i < len; i++) {
        sv->data[i] = IGRAPH_CALLOC(1, char);
        if (sv->data[i] == 0) {
            igraph_strvector_destroy(sv);
            IGRAPH_ERROR("strvector init failed", IGRAPH_ENOMEM);
        }
        sv->data[i][0] = '\0';
    }
    sv->len = len;

    return 0;
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
    long int i;
    IGRAPH_ASSERT(sv != 0);
    if (sv->data != 0) {
        for (i = 0; i < sv->len; i++) {
            if (sv->data[i] != 0) {
                IGRAPH_FREE(sv->data[i]);
            }
        }
        IGRAPH_FREE(sv->data);
    }
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
 * \param Pointer to a <type>char*</type>, the address of the string
 *   is stored here.
 *
 * Time complexity: O(1).
 */

void igraph_strvector_get(const igraph_strvector_t *sv, long int idx,
                          char **value) {
    IGRAPH_ASSERT(sv != 0);
    IGRAPH_ASSERT(sv->data != 0);
    IGRAPH_ASSERT(sv->data[idx] != 0);
    *value = sv->data[idx];
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

int igraph_strvector_set(igraph_strvector_t *sv, long int idx,
                         const char *value) {
    size_t value_len;

    IGRAPH_ASSERT(sv != 0);
    IGRAPH_ASSERT(sv->data != 0);

    value_len = strlen(value);
    if (sv->data[idx] == 0) {        
        sv->data[idx] = IGRAPH_CALLOC(value_len + 1, char);
        if (sv->data[idx] == 0) {
            IGRAPH_ERROR("strvector set failed", IGRAPH_ENOMEM);
        }
    } else {
        char *tmp = IGRAPH_REALLOC(sv->data[idx], value_len + 1, char);
        if (tmp == 0) {
            IGRAPH_ERROR("strvector set failed", IGRAPH_ENOMEM);
        }
        sv->data[idx] = tmp;
    }
    strcpy(sv->data[idx], value);

    return 0;
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
int igraph_strvector_set2(igraph_strvector_t *sv, long int idx,
                          const char *value, int len) {
    if (idx < 0 || idx >= sv->len) {
        IGRAPH_ERROR("String vector index out of bounds.", IGRAPH_EINVAL);
    }
    IGRAPH_ASSERT(sv != 0);
    IGRAPH_ASSERT(sv->data != 0);
    if (sv->data[idx] == 0) {
        sv->data[idx] = IGRAPH_CALLOC(len + 1, char);
        if (sv->data[idx] == 0) {
            IGRAPH_ERROR("strvector set failed", IGRAPH_ENOMEM);
        }
    } else {
        char *tmp = IGRAPH_REALLOC(sv->data[idx], (size_t) len + 1, char);
        if (tmp == 0) {
            IGRAPH_ERROR("strvector set failed", IGRAPH_ENOMEM);
        }
        sv->data[idx] = tmp;
    }
    memcpy(sv->data[idx], value, (size_t) len * sizeof(char));
    sv->data[idx][len] = '\0';

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_remove_section
 * \brief Removes a section from a string vector.
 * \todo repair realloc
 */

void igraph_strvector_remove_section(igraph_strvector_t *v, long int from,
                                     long int to) {
    long int i;
    /*   char **tmp; */

    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);

    for (i = from; i < to; i++) {
        if (v->data[i] != 0) {
            IGRAPH_FREE(v->data[i]);
        }
    }
    for (i = 0; i < v->len - to; i++) {
        v->data[from + i] = v->data[to + i];
    }

    v->len -= (to - from);

    /* try to make it smaller */
    /*   tmp=IGRAPH_REALLOC(v->data, v->len, char*); */
    /*   if (tmp!=0) { */
    /*     v->data=tmp; */
    /*   } */
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

void igraph_strvector_remove(igraph_strvector_t *v, long int elem) {
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);
    igraph_strvector_remove_section(v, elem, elem + 1);
}

/**
 * \ingroup strvector
 * \function igraph_strvector_move_interval
 * \brief Copies an interval of a string vector.
 */

void igraph_strvector_move_interval(igraph_strvector_t *v, long int begin,
                                    long int end, long int to) {
    long int i;
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);
    for (i = to; i < to + end - begin; i++) {
        if (v->data[i] != 0) {
            IGRAPH_FREE(v->data[i]);
        }
    }
    for (i = 0; i < end - begin; i++) {
        if (v->data[begin + i] != 0) {
            size_t len = strlen(v->data[begin + i]) + 1;
            v->data[to + i] = IGRAPH_CALLOC(len, char);
            memcpy(v->data[to + i], v->data[begin + i], sizeof(char)*len);
        }
    }
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

int igraph_strvector_copy(igraph_strvector_t *to,
                          const igraph_strvector_t *from) {
    long int i;
    char *str;
    IGRAPH_ASSERT(from != 0);
    /*   IGRAPH_ASSERT(from->data != 0); */
    to->data = IGRAPH_CALLOC(from->len, char*);
    if (to->data == 0) {
        IGRAPH_ERROR("Cannot copy string vector", IGRAPH_ENOMEM);
    }
    to->len = from->len;

    for (i = 0; i < from->len; i++) {
        int ret;
        igraph_strvector_get(from, i, &str);
        ret = igraph_strvector_set(to, i, str);
        if (ret != 0) {
            igraph_strvector_destroy(to);
            IGRAPH_ERROR("cannot copy string vector", ret);
        }
    }

    return 0;
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

int igraph_strvector_append(igraph_strvector_t *to,
                            const igraph_strvector_t *from) {
    long int len1 = igraph_strvector_size(to), len2 = igraph_strvector_size(from);
    long int i;
    igraph_bool_t error = 0;
    IGRAPH_CHECK(igraph_strvector_resize(to, len1 + len2));
    for (i = 0; i < len2; i++) {
        if (from->data[i][0] != '\0') {
            IGRAPH_FREE(to->data[len1 + i]);
            to->data[len1 + i] = strdup(from->data[i]);
            if (!to->data[len1 + i]) {
                error = 1;
                break;
            }
        }
    }
    if (error) {
        igraph_strvector_resize(to, len1);
        IGRAPH_ERROR("Cannot append string vector", IGRAPH_ENOMEM);
    }
    return 0;
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
    long int i, n = igraph_strvector_size(sv);
    char **tmp;

    for (i = 0; i < n; i++) {
        IGRAPH_FREE(sv->data[i]);
    }
    sv->len = 0;
    /* try to give back some memory */
    tmp = IGRAPH_REALLOC(sv->data, 1, char*);
    if (tmp != 0) {
        sv->data = tmp;
    }
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

int igraph_strvector_resize(igraph_strvector_t* v, long int newsize) {
    long int toadd = newsize - v->len, i, j;
    char **tmp;
    long int reallocsize = newsize;

    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);
    /*   printf("resize %li to %li\n", v->len, newsize); */
    if (newsize < v->len) {
        for (i = newsize; i < v->len; i++) {
            IGRAPH_FREE(v->data[i]);
        }
        /* try to give back some space */
        tmp = IGRAPH_REALLOC(v->data, (size_t) reallocsize, char*);
        /*     printf("resize %li to %li, %p\n", v->len, newsize, tmp); */
        if (tmp != 0) {
            v->data = tmp;
        }
    } else if (newsize > v->len) {
        igraph_bool_t error = 0;
        tmp = IGRAPH_REALLOC(v->data, (size_t) reallocsize, char*);
        if (tmp == 0) {
            IGRAPH_ERROR("cannot resize string vector", IGRAPH_ENOMEM);
        }
        v->data = tmp;

        for (i = 0; i < toadd; i++) {
            v->data[v->len + i] = IGRAPH_CALLOC(1, char);
            if (v->data[v->len + i] == 0) {
                error = 1;
                break;
            }
            v->data[v->len + i][0] = '\0';
        }
        if (error) {
            /* There was an error, free everything we've allocated so far */
            for (j = 0; j < i; j++) {
                if (v->data[v->len + i] != 0) {
                    IGRAPH_FREE(v->data[v->len + i]);
                }
            }
            /* Try to give back space */
            tmp = IGRAPH_REALLOC(v->data, (size_t) (v->len), char*);
            if (tmp != 0) {
                v->data = tmp;
            }
            IGRAPH_ERROR("Cannot resize string vector", IGRAPH_ENOMEM);
        }
    }
    v->len = newsize;

    return 0;
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

long int igraph_strvector_size(const igraph_strvector_t *sv) {
    IGRAPH_ASSERT(sv != 0);
    IGRAPH_ASSERT(sv->data != 0);
    return sv->len;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_add
 * \brief Adds an element to the back of a string vector.
 *
 * \param v The string vector.
 * \param value The string to add, it will be copied.
 * \return Error code.
 *
 * Time complexity: O(n+l), n is the total number of strings, l is the
 * length of the new string.
 */

int igraph_strvector_add(igraph_strvector_t *v, const char *value) {
    long int s = igraph_strvector_size(v);
    long int value_len = strlen(value);
    char **tmp;
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);
    tmp = IGRAPH_REALLOC(v->data, (size_t) s + 1, char*);
    if (tmp == 0) {
        IGRAPH_ERROR("cannot add string to string vector", IGRAPH_ENOMEM);
    }
    v->data = tmp;
    v->data[s] = IGRAPH_CALLOC(value_len + 1, char);
    if (v->data[s] == 0) {
        IGRAPH_ERROR("cannot add string to string vector", IGRAPH_ENOMEM);
    }
    strcpy(v->data[s], value);
    v->len += 1;

    return 0;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_permdelete
 * \brief Removes elements from a string vector (for internal use)
 */

void igraph_strvector_permdelete(igraph_strvector_t *v, const igraph_vector_t *index,
                                 long int nremove) {
    long int i;
    char **tmp;
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);

    for (i = 0; i < igraph_strvector_size(v); i++) {
        if (VECTOR(*index)[i] != 0) {
            v->data[ (long int) VECTOR(*index)[i] - 1 ] = v->data[i];
        } else {
            IGRAPH_FREE(v->data[i]);
        }
    }
    /* Try to make it shorter */
    tmp = IGRAPH_REALLOC(v->data, v->len - nremove ?
                         (size_t) (v->len - nremove) : 1, char*);
    if (tmp != 0) {
        v->data = tmp;
    }
    v->len -= nremove;
}

/**
 * \ingroup strvector
 * \function igraph_strvector_remove_negidx
 * \brief Removes elements from a string vector (for internal use)
 */

void igraph_strvector_remove_negidx(igraph_strvector_t *v, const igraph_vector_t *neg,
                                    long int nremove) {
    long int i, idx = 0;
    char **tmp;
    IGRAPH_ASSERT(v != 0);
    IGRAPH_ASSERT(v->data != 0);
    for (i = 0; i < igraph_strvector_size(v); i++) {
        if (VECTOR(*neg)[i] >= 0) {
            v->data[idx++] = v->data[i];
        } else {
            IGRAPH_FREE(v->data[i]);
        }
    }
    /* Try to give back some memory */
    tmp = IGRAPH_REALLOC(v->data, v->len - nremove ?
                         (size_t) (v->len - nremove) : 1, char*);
    if (tmp != 0) {
        v->data = tmp;
    }
    v->len -= nremove;
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
int igraph_strvector_print(const igraph_strvector_t *v, FILE *file,
                           const char *sep) {

    long int i, n = igraph_strvector_size(v);
    if (n != 0) {
        fprintf(file, "%s", STR(*v, 0));
    }
    for (i = 1; i < n; i++) {
        fprintf(file, "%s%s", sep, STR(*v, i));
    }
    return IGRAPH_SUCCESS;
}

int igraph_strvector_index(const igraph_strvector_t *v,
                           igraph_strvector_t *newv,
                           const igraph_vector_t *idx) {

    long int i, newlen = igraph_vector_size(idx);
    IGRAPH_CHECK(igraph_strvector_resize(newv, newlen));

    for (i = 0; i < newlen; i++) {
        long int j = (long int) VECTOR(*idx)[i];
        char *str;
        igraph_strvector_get(v, j, &str);
        igraph_strvector_set(newv, i, str);
    }

    return 0;
}
