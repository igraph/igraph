/* -*- mode: C -*-  */
/*
   IGraph library.
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
#include "igraph_types_internal.h"

int main() {

    igraph_trie_t trie;
    long int id;
    int i;
    char *str;

    /* init */
    igraph_trie_init(&trie, 0);

    /* add and get values */
    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %li\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %li\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %li\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %li\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %li\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %li\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %li\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %li\n", id);

    igraph_trie_get(&trie, "a", &id);
    printf("a:     %li\n", id);
    igraph_trie_get(&trie, "axon", &id);
    printf("axon:  %li\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %li\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %li\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %li\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %li\n", id);

    /* check for existence */
    igraph_trie_check(&trie, "head", &id);
    printf("head:  %li\n", id);
    igraph_trie_check(&trie, "alma", &id);
    printf("alma:  %li\n", id);

    /* destroy */
    igraph_trie_destroy(&trie);

    /* the same with index */
    igraph_trie_init(&trie, 1);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %li\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %li\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %li\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %li\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %li\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %li\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %li\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %li\n", id);

    igraph_trie_get(&trie, "a", &id);
    printf("a:     %li\n", id);
    igraph_trie_get(&trie, "axon", &id);
    printf("axon:  %li\n", id);

    igraph_trie_get(&trie, "hello", &id);
    printf("hello: %li\n", id);
    igraph_trie_get(&trie, "hepp", &id);
    printf("hepp:  %li\n", id);
    igraph_trie_get(&trie, "alma", &id);
    printf("alma:  %li\n", id);
    igraph_trie_get(&trie, "also", &id);
    printf("also:  %li\n", id);

    /* check for existence */
    igraph_trie_check(&trie, "head", &id);
    printf("head:  %li\n", id);
    igraph_trie_check(&trie, "alma", &id);
    printf("alma:  %li\n", id);

    for (i = 0; i < igraph_trie_size(&trie); i++) {
        igraph_trie_idx(&trie, i, &str);
        printf("%d: %s\n", i, str);
    }
    igraph_trie_destroy(&trie);

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        return 1;
    }

    return 0;
}
