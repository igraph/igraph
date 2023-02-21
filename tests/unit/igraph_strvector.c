/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>

#include "test_utilities.h"

void strvector_print(const igraph_strvector_t *sv) {
    igraph_integer_t i, s = igraph_strvector_size(sv);
    for (i = 0; i < s; i++) {
        printf("\"%s\"\n", STR(*sv, i));
    }
    printf("\n");
}

int main(void) {

    igraph_strvector_t sv1, sv2, sv3, sv4;

    printf("igraph_strvector_init, igraph_strvector_destroy\n");
    igraph_strvector_init(&sv1, 10);
    igraph_strvector_destroy(&sv1);
    igraph_strvector_init(&sv1, 0);
    igraph_strvector_destroy(&sv1);

    printf("igraph_strvector_get, igraph_strvector_set\n");
    igraph_strvector_init(&sv1, 5);
    strvector_print(&sv1);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    strvector_print(&sv1);

    printf("igraph_strvector_remove_section\n");
    igraph_strvector_remove_section(&sv1, 0, 5);
    if (igraph_strvector_size(&sv1) != 0) {
        return 1;
    }
    printf("resize to 10 and then back to 5\n");
    igraph_strvector_resize(&sv1, 10);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    igraph_strvector_resize(&sv1, 5);
    strvector_print(&sv1);
    printf("resize to 0\n");
    igraph_strvector_resize(&sv1, 0);
    if (igraph_strvector_size(&sv1) != 0) {
        return 1;
    }
    printf("resize to 10 and then back to 5 again\n");
    igraph_strvector_resize(&sv1, 10);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");
    igraph_strvector_set(&sv1, 3, "three");
    igraph_strvector_set(&sv1, 4, "four");
    igraph_strvector_resize(&sv1, 5);
    strvector_print(&sv1);

    printf("igraph_strvector_copy\n");
    igraph_strvector_init_copy(&sv2, &sv1);
    strvector_print(&sv1);

    igraph_strvector_resize(&sv1, 0);
    igraph_strvector_destroy(&sv2);
    igraph_strvector_init_copy(&sv2, &sv1);
    if (igraph_strvector_size(&sv2) != 0) {
        return 2;
    }
    igraph_strvector_destroy(&sv2);

    printf("igraph_strvector_push_back\n");
    igraph_strvector_push_back(&sv1, "zeroth");
    igraph_strvector_push_back(&sv1, "first");
    igraph_strvector_push_back(&sv1, "second");
    igraph_strvector_push_back(&sv1, "third");
    igraph_strvector_push_back(&sv1, "fourth");
    strvector_print(&sv1);

    printf("igraph_strvector_push_back_len\n");
    igraph_strvector_push_back_len(&sv1, "extra zeroth", 5);
    igraph_strvector_push_back_len(&sv1, "extra first", 100);
    igraph_strvector_push_back_len(&sv1, "extra second", 0);
    strvector_print(&sv1);
    igraph_strvector_destroy(&sv1);

    printf("igraph_strvector_append\n");
    printf("===\n");
    igraph_strvector_init(&sv1, 0);
    igraph_strvector_init(&sv2, 0);
    igraph_strvector_append(&sv1, &sv2);
    strvector_print(&sv1);
    printf("===\n");

    igraph_strvector_resize(&sv1, 3);
    igraph_strvector_append(&sv1, &sv2);
    strvector_print(&sv1);
    printf("===\n");

    igraph_strvector_append(&sv2, &sv1);
    strvector_print(&sv2);
    printf("===\n");

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

    printf("igraph_strvector_merge\n");
    igraph_strvector_init(&sv1, 3);
    igraph_strvector_set(&sv1, 0, "zero");
    igraph_strvector_set(&sv1, 1, "one");
    igraph_strvector_set(&sv1, 2, "two");

    igraph_strvector_init(&sv2, 2);
    igraph_strvector_set(&sv2, 0, "a");
    igraph_strvector_set(&sv2, 1, "b");

    igraph_strvector_init_copy(&sv3, &sv1);
    igraph_strvector_init_copy(&sv4, &sv2);

    igraph_strvector_merge(&sv1, &sv2);
    IGRAPH_ASSERT(igraph_strvector_size(&sv2) == 0);

    igraph_strvector_append(&sv3, &sv4);
    IGRAPH_ASSERT(igraph_strvector_size(&sv1) == igraph_strvector_size(&sv3));

    for (igraph_integer_t i=0; i < igraph_strvector_size(&sv1); ++i) {
        IGRAPH_ASSERT(strcmp(STR(sv1, i), STR(sv3, i)) == 0);
    }

    igraph_strvector_destroy(&sv1);
    igraph_strvector_destroy(&sv2);
    igraph_strvector_destroy(&sv3);
    igraph_strvector_destroy(&sv4);

    printf("clear\n");
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

    VERIFY_FINALLY_STACK();
    return 0;
}
