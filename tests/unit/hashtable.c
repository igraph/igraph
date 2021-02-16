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

#include <igraph.h>
#include "core/hashtable.h"

#include "test_utilities.inc"

int main() {

    igraph_hashtable_t ht;
    char *str;
    const igraph_strvector_t *keys;
    long int i;

    /* init and destroy */
    igraph_hashtable_init(&ht);
    igraph_hashtable_destroy(&ht);

    /* init, add some elements and destroy */
    igraph_hashtable_init(&ht);
    igraph_hashtable_addset(&ht, "color", "green", "red");
    igraph_hashtable_addset(&ht, "size", "", "4");
    igraph_hashtable_addset(&ht, "color", "", "grey");
    igraph_hashtable_addset(&ht, "shape", "", "circle");
    igraph_hashtable_addset(&ht, "shape", "", "diamond");
    igraph_hashtable_destroy(&ht);

    /* reset */
    igraph_hashtable_init(&ht);
    igraph_hashtable_addset(&ht, "color", "green", "red");
    igraph_hashtable_addset(&ht, "size", "", "4");
    igraph_hashtable_addset(&ht, "color", "", "grey");
    igraph_hashtable_addset(&ht, "shape", "", "circle");
    igraph_hashtable_addset(&ht, "shape", "", "diamond");
    igraph_hashtable_reset(&ht);
    igraph_hashtable_addset(&ht, "color", "green", "red");
    igraph_hashtable_addset(&ht, "size", "", "4");
    igraph_hashtable_addset(&ht, "color", "", "grey");
    igraph_hashtable_addset(&ht, "shape", "", "circle");
    igraph_hashtable_addset(&ht, "shape", "", "diamond");
    igraph_hashtable_destroy(&ht);

    /* Check semantics */
    igraph_hashtable_init(&ht);
    igraph_hashtable_addset(&ht, "color", "green", "red");
    igraph_hashtable_addset(&ht, "size", "", "4");
    igraph_hashtable_addset(&ht, "color", "", "grey");
    igraph_hashtable_addset(&ht, "shape", "", "circle");
    igraph_hashtable_addset(&ht, "shape", "", "diamond");

    igraph_hashtable_get(&ht, "color", &str);
    printf("color: %s\n", str);
    igraph_hashtable_get(&ht, "size", &str);
    printf("size: %s\n", str);
    igraph_hashtable_get(&ht, "shape", &str);
    printf("shape: %s\n", str);

    igraph_hashtable_reset(&ht);

    igraph_hashtable_get(&ht, "color", &str);
    printf("color: %s\n", str);
    igraph_hashtable_get(&ht, "size", &str);
    printf("size: %s\n", str);
    igraph_hashtable_get(&ht, "shape", &str);
    printf("shape: %s\n", str);

    igraph_hashtable_getkeys(&ht, &keys);
    for (i = 0; i < igraph_strvector_size(keys); i++) {
        igraph_strvector_get(keys, i, &str);
        printf("%s ", str);
    }
    printf("\n");

    igraph_hashtable_destroy(&ht);

    /* addset2 */
    igraph_hashtable_init(&ht);
    igraph_hashtable_addset2(&ht, "color", "green", "redddd", 3);
    igraph_hashtable_addset2(&ht, "size", "", "4111", 1);
    igraph_hashtable_addset2(&ht, "color", "", "greysdsdf", 4);
    igraph_hashtable_addset2(&ht, "shape", "", "circle", 6);
    igraph_hashtable_addset(&ht, "shape", "", "diamond");

    igraph_hashtable_get(&ht, "color", &str);
    printf("color: %s\n", str);
    igraph_hashtable_get(&ht, "size", &str);
    printf("size: %s\n", str);
    igraph_hashtable_get(&ht, "shape", &str);
    printf("shape: %s\n", str);

    igraph_hashtable_reset(&ht);

    igraph_hashtable_get(&ht, "color", &str);
    printf("color: %s\n", str);
    igraph_hashtable_get(&ht, "size", &str);
    printf("size: %s\n", str);
    igraph_hashtable_get(&ht, "shape", &str);
    printf("shape: %s\n", str);

    igraph_hashtable_getkeys(&ht, &keys);
    for (i = 0; i < igraph_strvector_size(keys); i++) {
        igraph_strvector_get(keys, i, &str);
        printf("%s ", str);
    }
    printf("\n");

    igraph_hashtable_destroy(&ht);

    VERIFY_FINALLY_STACK();

    return 0;
}
