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

#include <igraph.h>

void strvector_print(const igraph_strvector_t sv) {
    igraph_int_t i, s = igraph_strvector_size(&sv);
    for (i = 0; i < s; i++) {
        printf("---%s---\n", igraph_strvector_get(&sv, i));
    }
    printf("\n");
}

int main(void) {

    igraph_strvector_t sv1, sv2;

    /* Initialize the library. */
    igraph_setup();

    printf("Initializing and setting elements:\n");
    igraph_strvector_init(&sv1, 5);

    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    strvector_print(sv1);

    printf("strvector size after first removing one element, and then the rest:\n");
    igraph_strvector_remove(&sv1, 4);
    igraph_strvector_remove_section(&sv1, 0, 4);
    printf("%" IGRAPH_PRId "\n\n", igraph_strvector_size(&sv1));

    printf("Resize to three elements, and set them:\n");
    igraph_strvector_resize(&sv1, 3);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    strvector_print(sv1);

    printf("Then copy the first element over the third element:\n");
    igraph_strvector_set(&sv1, 2, igraph_strvector_get(&sv1, 0));
    strvector_print(sv1);

    printf("Make a copy of the strvector and set the last element of the copy:\n");
    igraph_strvector_init_copy(&sv2, &sv1);
    igraph_strvector_set(&sv2, 2, "copy two");
    strvector_print(sv2);

    printf("Append the copy to the strvector:\n");
    igraph_strvector_append(&sv1, &sv2);
    strvector_print(sv1);

    printf("Add two strings at the end:\n");
    igraph_strvector_push_back(&sv1, "zeroth");
    igraph_strvector_push_back(&sv1, "first");
    strvector_print(sv1);
    igraph_strvector_destroy(&sv1);

    printf("strvector size after clearing it:\n");
    igraph_strvector_clear(&sv2);
    printf("%" IGRAPH_PRId "\n", igraph_strvector_size(&sv2));

    igraph_strvector_destroy(&sv2);

    return 0;
}
