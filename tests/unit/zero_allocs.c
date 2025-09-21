/*
   igraph library.
   Copyright (C) 2021-2024  The igraph development team <igraph@igraph.org>

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
    int *a;
    char *b;

    /* Macros */

    a = IGRAPH_CALLOC(0, int);
    IGRAPH_ASSERT(a);

    a = IGRAPH_REALLOC(a, 0, int);
    IGRAPH_ASSERT(a);

    IGRAPH_FREE(a);
    IGRAPH_ASSERT(!a); /* IGRAPH_FREE(a) sets 'a' to NULL */

    /* We use type 'char' to work around warnings from some moderns compilers such as Clang 22:
     * error: allocation of insufficient size '1' for type 'int' with size '4' [-Werror,-Walloc-size] */
    b = IGRAPH_MALLOC(0);
    IGRAPH_ASSERT(b);

    IGRAPH_FREE(b);
    IGRAPH_ASSERT(!b);

    /* Functions */

    a = igraph_calloc(0, sizeof(*a));
    IGRAPH_ASSERT(a);

    a = igraph_realloc(a, 0);
    IGRAPH_ASSERT(a);

    igraph_free(a);

    a = igraph_malloc(0);
    IGRAPH_ASSERT(a);

    igraph_free(a);

    VERIFY_FINALLY_STACK();
    return 0;
}
