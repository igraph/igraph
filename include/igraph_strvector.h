/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_STRVECTOR_H
#define IGRAPH_STRVECTOR_H

#include "igraph_decls.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/**
 * Vector of strings
 * \ingroup internal
 */

typedef struct s_igraph_strvector {
    char **data;
    long int len;
} igraph_strvector_t;

/**
 * \define STR
 * Indexing string vectors
 *
 * This is a macro which allows to query the elements of a string vector in
 * simpler way than \ref igraph_strvector_get(). Note this macro cannot be
 * used to set an element, for that use \ref igraph_strvector_set().
 * \param sv The string vector
 * \param i The the index of the element.
 * \return The element at position \p i.
 *
 * Time complexity: O(1).
 */
#define STR(sv,i) ((const char *)((sv).data[(i)]))

#define IGRAPH_STRVECTOR_NULL { 0,0 }
#define IGRAPH_STRVECTOR_INIT_FINALLY(v, size) \
    do { IGRAPH_CHECK(igraph_strvector_init(v, size)); \
        IGRAPH_FINALLY( (igraph_finally_func_t*) igraph_strvector_destroy, v); } while (0)

DECLDIR int igraph_strvector_init(igraph_strvector_t *sv, long int len);
DECLDIR void igraph_strvector_destroy(igraph_strvector_t *sv);
DECLDIR long int igraph_strvector_size(const igraph_strvector_t *sv);
DECLDIR void igraph_strvector_get(const igraph_strvector_t *sv,
                                  long int idx, char **value);
DECLDIR int igraph_strvector_set(igraph_strvector_t *sv, long int idx,
                                 const char *value);
DECLDIR int igraph_strvector_set2(igraph_strvector_t *sv, long int idx,
                                  const char *value, int len);
DECLDIR void igraph_strvector_clear(igraph_strvector_t *sv);
DECLDIR void igraph_strvector_remove_section(igraph_strvector_t *v, long int from,
        long int to);
DECLDIR void igraph_strvector_remove(igraph_strvector_t *v, long int elem);
DECLDIR void igraph_strvector_move_interval(igraph_strvector_t *v, long int begin,
        long int end, long int to);
DECLDIR int igraph_strvector_copy(igraph_strvector_t *to,
                                  const igraph_strvector_t *from);
DECLDIR int igraph_strvector_append(igraph_strvector_t *to,
                                    const igraph_strvector_t *from);
DECLDIR int igraph_strvector_resize(igraph_strvector_t* v, long int newsize);
DECLDIR int igraph_strvector_add(igraph_strvector_t *v, const char *value);
DECLDIR void igraph_strvector_permdelete(igraph_strvector_t *v, const igraph_vector_t *index,
        long int nremove);
DECLDIR void igraph_strvector_remove_negidx(igraph_strvector_t *v, const igraph_vector_t *neg,
        long int nremove);
DECLDIR int igraph_strvector_print(const igraph_strvector_t *v, FILE *file,
                                   const char *sep);

DECLDIR int igraph_strvector_index(const igraph_strvector_t *v,
                                   igraph_strvector_t *newv,
                                   const igraph_vector_t *idx);


__END_DECLS

#endif
