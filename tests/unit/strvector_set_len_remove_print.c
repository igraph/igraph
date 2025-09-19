/*
   igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_strvector_t sv;
    char *test_string = "This is a string.";
    char *test_string2 = "A completely different one.";

    printf("Two strings in a vector.\n");
    igraph_strvector_init(&sv, 5);
    IGRAPH_ASSERT(igraph_strvector_set_len(&sv, 0, test_string, strlen(test_string)) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_strvector_set_len(&sv, 4, test_string2, strlen(test_string2)) == IGRAPH_SUCCESS);
    igraph_strvector_print(&sv, " | ");

    printf("\nRemove a nonexistent one.\n");
    igraph_strvector_remove(&sv, 1);
    igraph_strvector_print(&sv, " | ");

    printf("\nRemove one.\n");
    igraph_strvector_remove(&sv, 0);
    igraph_strvector_print(&sv, " | ");
    igraph_strvector_destroy(&sv);

    printf("\nOverwriting a string.\n");
    igraph_strvector_init(&sv, 5);
    IGRAPH_ASSERT(igraph_strvector_set_len(&sv, 2, test_string2, strlen(test_string2)) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_strvector_set_len(&sv, 2, test_string, strlen(test_string)) == IGRAPH_SUCCESS);
    igraph_strvector_print(&sv, " | ");
    igraph_strvector_destroy(&sv);

    printf("\n");

    VERIFY_FINALLY_STACK();
    return 0;
}
