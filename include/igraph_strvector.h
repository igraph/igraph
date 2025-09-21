/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_STRVECTOR_H
#define IGRAPH_STRVECTOR_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_vector.h"

IGRAPH_BEGIN_C_DECLS

/**
 * Vector of strings
 * \ingroup internal
 */

typedef struct s_igraph_strvector {
    /* Empty strings "" are represented using NULL. */
    const char **stor_begin;
    const char **stor_end;
    const char **end;
} igraph_strvector_t;

/**
 * \define STR
 * \brief Indexing string vectors.
 *
 * This is a macro that allows to query the elements of a string vector, just
 * like \ref igraph_strvector_get(). Note this macro cannot be used to set an
 * element. Use \ref igraph_strvector_set() to set an element instead.
 *
 * \param sv The string vector
 * \param i The index of the element.
 * \return The element at position \p i.
 *
 * Time complexity: O(1).
 *
 * \deprecated-by igraph_strvector_get 0.10.9
 */
#define STR(sv,i) \
    (IGRAPH_PREPROCESSOR_WARNING("STR() is deprecated. Use igraph_strvector_get() instead.") \
     igraph_strvector_get(&sv, i))

#define IGRAPH_STRVECTOR_NULL { 0,0,0 }
#define IGRAPH_STRVECTOR_INIT_FINALLY(sv, size) \
    do { IGRAPH_CHECK(igraph_strvector_init(sv, size)); \
        IGRAPH_FINALLY( igraph_strvector_destroy, sv); } while (0)

IGRAPH_EXPORT igraph_error_t igraph_strvector_init(igraph_strvector_t *sv, igraph_int_t len);
IGRAPH_EXPORT void igraph_strvector_destroy(igraph_strvector_t *sv);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_strvector_size(const igraph_strvector_t *sv);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t igraph_strvector_capacity(const igraph_strvector_t *sv);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE const char *igraph_strvector_get(const igraph_strvector_t *sv, igraph_int_t idx);
IGRAPH_EXPORT igraph_error_t igraph_strvector_set(
    igraph_strvector_t *sv, igraph_int_t idx, const char *value);
IGRAPH_EXPORT igraph_error_t igraph_strvector_set_len(
    igraph_strvector_t *sv, igraph_int_t idx, const char *value, size_t len);
IGRAPH_EXPORT void igraph_strvector_clear(igraph_strvector_t *sv);
IGRAPH_EXPORT void igraph_strvector_remove_section(
    igraph_strvector_t *v, igraph_int_t from, igraph_int_t to);
IGRAPH_EXPORT void igraph_strvector_remove(
    igraph_strvector_t *v, igraph_int_t elem);
IGRAPH_EXPORT igraph_error_t igraph_strvector_init_copy(
    igraph_strvector_t *to, const igraph_strvector_t *from);
IGRAPH_EXPORT igraph_error_t igraph_strvector_append(
    igraph_strvector_t *to, const igraph_strvector_t *from);
IGRAPH_EXPORT igraph_error_t igraph_strvector_merge(
    igraph_strvector_t *to, igraph_strvector_t *from);
IGRAPH_EXPORT void igraph_strvector_swap(igraph_strvector_t *v1, igraph_strvector_t *v2);
IGRAPH_EXPORT igraph_error_t igraph_strvector_update(
    igraph_strvector_t *to, const igraph_strvector_t *from);
IGRAPH_EXPORT igraph_error_t igraph_strvector_resize(
    igraph_strvector_t* v, igraph_int_t newsize);
IGRAPH_EXPORT void igraph_strvector_resize_min(igraph_strvector_t *sv);
IGRAPH_EXPORT igraph_error_t igraph_strvector_push_back(igraph_strvector_t *v,
        const char *value);
IGRAPH_EXPORT igraph_error_t igraph_strvector_push_back_len(igraph_strvector_t *v,
        const char *value, size_t len);
IGRAPH_EXPORT igraph_error_t igraph_strvector_fprint(const igraph_strvector_t *v, FILE *file,
                                                     const char *sep);
IGRAPH_EXPORT igraph_error_t igraph_strvector_print(const igraph_strvector_t *v,
                                                    const char *sep);

IGRAPH_EXPORT igraph_error_t igraph_strvector_index(const igraph_strvector_t *v,
                                         igraph_strvector_t *newv,
                                         const igraph_vector_int_t *idx);

IGRAPH_EXPORT igraph_error_t igraph_strvector_reserve(igraph_strvector_t *sv,
                                                      igraph_int_t capacity);

IGRAPH_EXPORT void igraph_strvector_swap_elements(igraph_strvector_t *sv,
                                                  igraph_int_t i, igraph_int_t j);

IGRAPH_END_C_DECLS

#endif
