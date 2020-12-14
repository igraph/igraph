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

#ifndef IGRAPH_CORE_HASHTABLE_H
#define IGRAPH_CORE_HASHTABLE_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
    #define __BEGIN_DECLS extern "C" {
    #define __END_DECLS }
#else
    #define __BEGIN_DECLS /* empty */
    #define __END_DECLS /* empty */
#endif

#include "igraph_types.h"
#include "igraph_strvector.h"

#include "core/trie.h"

__BEGIN_DECLS

/* string -> string hash table */

typedef struct igraph_hashtable_t {
    igraph_trie_t keys;
    igraph_strvector_t elements;
    igraph_strvector_t defaults;
} igraph_hashtable_t;

int igraph_hashtable_init(igraph_hashtable_t *ht);
void igraph_hashtable_destroy(igraph_hashtable_t *ht);
int igraph_hashtable_addset(igraph_hashtable_t *ht,
                            const char *key, const char *def,
                            const char *elem);
int igraph_hashtable_addset2(igraph_hashtable_t *ht,
                             const char *key, const char *def,
                             const char *elem, int elemlen);
int igraph_hashtable_get(igraph_hashtable_t *ht,
                         const char *key, char **elem);
int igraph_hashtable_getkeys(igraph_hashtable_t *ht,
                             const igraph_strvector_t **sv);
int igraph_hashtable_reset(igraph_hashtable_t *ht);

__END_DECLS

#endif
