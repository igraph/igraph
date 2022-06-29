/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2020  The igraph development team

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

#ifndef IGRAPH_CORE_TRIE_H
#define IGRAPH_CORE_TRIE_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_strvector.h"
#include "igraph_vector.h"
#include "igraph_vector_ptr.h"

__BEGIN_DECLS

/**
 * Trie data type
 * \ingroup internal
 */

typedef struct s_igraph_trie_node {
    igraph_strvector_t strs;
    igraph_vector_ptr_t children;
    igraph_vector_int_t values;
} igraph_trie_node_t;

typedef struct s_igraph_trie {
    igraph_trie_node_t node;
    igraph_integer_t maxvalue;
    igraph_bool_t storekeys;
    igraph_strvector_t keys;
} igraph_trie_t;

#define IGRAPH_TRIE_NULL \
        { { IGRAPH_STRVECTOR_NULL, IGRAPH_VECTOR_PTR_NULL, IGRAPH_VECTOR_NULL}, \
            0, 0, IGRAPH_STRVECTOR_NULL }
#define IGRAPH_TRIE_INIT_FINALLY(tr, sk) \
    do { IGRAPH_CHECK(igraph_trie_init(tr, sk)); \
        IGRAPH_FINALLY(igraph_trie_destroy, tr); } while (0)

IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_trie_init(igraph_trie_t *t, igraph_bool_t storekeys);
IGRAPH_PRIVATE_EXPORT void igraph_trie_destroy(igraph_trie_t *t);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_trie_get(igraph_trie_t *t, const char *key, igraph_integer_t *id);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_trie_check(igraph_trie_t *t, const char *key, igraph_integer_t *id);
IGRAPH_PRIVATE_EXPORT igraph_error_t igraph_trie_get_len(igraph_trie_t *t, const char *key, igraph_integer_t length,
                                           igraph_integer_t *id);
IGRAPH_PRIVATE_EXPORT const char* igraph_trie_idx(igraph_trie_t *t, igraph_integer_t idx);
IGRAPH_PRIVATE_EXPORT igraph_integer_t igraph_trie_size(igraph_trie_t *t);

const igraph_strvector_t* igraph_i_trie_borrow_keys(igraph_trie_t *t);

__END_DECLS

#endif
