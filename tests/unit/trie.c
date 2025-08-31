/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

#include <stdio.h>
#include <igraph.h>

#include "core/trie.h"

#include "test_utilities.h"

int main(void) {

    igraph_trie_t trie;
    igraph_int_t id;
    igraph_int_t i;
    const char *str;

    /* init */
    igraph_trie_init(&trie, 0);

    /* add and get values */
    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %" IGRAPH_PRId "\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %" IGRAPH_PRId "\n", id);

    igraph_trie_get(&trie, "a", &id);
    printf("a:     %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "axon", &id);
    printf("axon:  %" IGRAPH_PRId "\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %" IGRAPH_PRId "\n", id);

    /* check for existence */
    igraph_trie_check(&trie, "head", &id);
    printf("head:  %" IGRAPH_PRId "\n", id);
    igraph_trie_check(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);

    /* destroy */
    igraph_trie_destroy(&trie);

    /* the same with index */
    igraph_trie_init(&trie, 1);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %" IGRAPH_PRId "\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %" IGRAPH_PRId "\n", id);

    igraph_trie_get(&trie, "a", &id);
    printf("a:     %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "axon", &id);
    printf("axon:  %" IGRAPH_PRId "\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %" IGRAPH_PRId "\n", id);

    /* check for existence */
    igraph_trie_check(&trie, "head", &id);
    printf("head:  %" IGRAPH_PRId "\n", id);
    igraph_trie_check(&trie, "alma", &id);
    printf("alma:  %" IGRAPH_PRId "\n", id);

    for (i = 0; i < igraph_trie_size(&trie); i++) {
        str = igraph_trie_idx(&trie, i);
        printf("%" IGRAPH_PRId ": %s\n", i, str);
    }

    /* prevent insertion of empty key */
    igraph_set_error_handler(igraph_error_handler_ignore);
    IGRAPH_ASSERT(igraph_trie_get(&trie, "", &id) == IGRAPH_EINVAL);

    igraph_trie_destroy(&trie);

    VERIFY_FINALLY_STACK();

    return 0;
}
