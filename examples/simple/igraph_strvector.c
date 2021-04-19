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

void strvector_print(const igraph_strvector_t *sv) {
    long int i, s = igraph_strvector_size(sv);
    for (i = 0; i < s; i++) {
        printf("---%s---\n", STR(*sv, i));
    }
}

int main() {

    igraph_strvector_t sv1, sv2;
    char *str1;
    int i;

    /* igraph_strvector_init, igraph_strvector_destroy */
    igraph_strvector_init(&sv1, 10);
    igraph_strvector_destroy(&sv1);
    igraph_strvector_init(&sv1, 0);
    igraph_strvector_destroy(&sv1);

    /* igraph_strvector_get, igraph_strvector_set */
    igraph_strvector_init(&sv1, 5);
    for (i = 0; i < igraph_strvector_size(&sv1); i++) {
        igraph_strvector_get(&sv1, i, &str1);
        printf("---%s---\n", str1);
    }
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    for (i = 0; i < igraph_strvector_size(&sv1); i++) {
        igraph_strvector_get(&sv1, i, &str1);
        printf("---%s---\n", str1);
    }

    /* igraph_strvector_remove_section, igraph_strvector_remove,
       igraph_strvector_resize, igraph_strvector_size */
    igraph_strvector_remove_section(&sv1, 0, 5);
    if (igraph_strvector_size(&sv1) != 0) {
        return 1;
    }
    igraph_strvector_resize(&sv1, 10);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    igraph_strvector_resize(&sv1, 5);
    for (i = 0; i < igraph_strvector_size(&sv1); i++) {
        igraph_strvector_get(&sv1, i, &str1);
        printf("---%s---\n", str1);
    }
    igraph_strvector_resize(&sv1, 0);
    if (igraph_strvector_size(&sv1) != 0) {
        return 1;
    }
    igraph_strvector_resize(&sv1, 10);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    igraph_strvector_resize(&sv1, 5);
    for (i = 0; i < igraph_strvector_size(&sv1); i++) {
        igraph_strvector_get(&sv1, i, &str1);
        printf("---%s---\n", str1);
    }

    /* igraph_strvector_move_interval */
    igraph_strvector_move_interval(&sv1, 3, 5, 0);
    for (i = 0; i < igraph_strvector_size(&sv1); i++) {
        igraph_strvector_get(&sv1, i, &str1);
        printf("---%s---\n", str1);
    }

    /* igraph_strvector_copy */
    igraph_strvector_copy(&sv2, &sv1);
    for (i = 0; i < igraph_strvector_size(&sv2); i++) {
        igraph_strvector_get(&sv2, i, &str1);
        printf("---%s---\n", str1);
    }
    igraph_strvector_resize(&sv1, 0);
    igraph_strvector_destroy(&sv2);
    igraph_strvector_copy(&sv2, &sv1);
    if (igraph_strvector_size(&sv2) != 0) {
        return 2;
    }
    igraph_strvector_destroy(&sv2);

    /* igraph_strvector_add */
    igraph_strvector_add(&sv1, "zeroth");
    igraph_strvector_add(&sv1, "first");
    igraph_strvector_add(&sv1, "second");
    igraph_strvector_add(&sv1, "third");
    igraph_strvector_add(&sv1, "fourth");
    for (i = 0; i < igraph_strvector_size(&sv1); i++) {
        igraph_strvector_get(&sv1, i, &str1);
        printf("---%s---\n", str1);
    }

    /* TODO: igraph_strvector_permdelete */
    /* TODO: igraph_strvector_remove_negidx */

    igraph_strvector_destroy(&sv1);

    /* append */
    printf("---\n");
    igraph_strvector_init(&sv1, 0);
    igraph_strvector_init(&sv2, 0);
    igraph_strvector_append(&sv1, &sv2);
    strvector_print(&sv1);
    printf("---\n");

    igraph_strvector_resize(&sv1, 3);
    igraph_strvector_append(&sv1, &sv2);
    strvector_print(&sv1);
    printf("---\n");

    igraph_strvector_append(&sv2, &sv1);
    strvector_print(&sv2);
    printf("---\n");

    igraph_strvector_set(&sv1, 0, "0");
    igraph_strvector_set(&sv1, 1, "1");
    igraph_strvector_set(&sv1, 2, "2");
    igraph_strvector_set(&sv2, 0, "3");
    igraph_strvector_set(&sv2, 1, "4");
    igraph_strvector_set(&sv2, 2, "5");
    igraph_strvector_append(&sv1, &sv2);
    strvector_print(&sv1);

    igraph_strvector_destroy(&sv1);
    igraph_strvector_destroy(&sv2);

    /* clear */
    igraph_strvector_init(&sv1, 3);
    igraph_strvector_set(&sv1, 0, "0");
    igraph_strvector_set(&sv1, 1, "1");
    igraph_strvector_set(&sv1, 2, "2");
    igraph_strvector_clear(&sv1);
    if (igraph_strvector_size(&sv1) != 0) {
        return 3;
    }
    igraph_strvector_resize(&sv1, 4);
    strvector_print(&sv1);
    igraph_strvector_set(&sv1, 0, "one");
    igraph_strvector_set(&sv1, 2, "two");
    strvector_print(&sv1);
    igraph_strvector_destroy(&sv1);

    /* STR */

    igraph_strvector_init(&sv1, 5);
    igraph_strvector_set(&sv1, 0, "one");
    igraph_strvector_set(&sv1, 1, "two");
    igraph_strvector_set(&sv1, 2, "three");
    igraph_strvector_set(&sv1, 3, "four");
    igraph_strvector_set(&sv1, 4, "five");
    strvector_print(&sv1);
    igraph_strvector_destroy(&sv1);

    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        return 4;
    }

    return 0;
}
