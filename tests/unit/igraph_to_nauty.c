/*
   IGraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

void print(igraph_t *g) {
    char *str = NULL;

    igraph_to_nauty(g, &str);
    printf("%s\n", str);


    IGRAPH_FREE(str);
    igraph_destroy(g);
}


int main(void) {
    igraph_t g;

    printf("Should be Gr`HOk:\n");
    igraph_small(&g, 8, IGRAPH_UNDIRECTED,
                 0,1, 0,2, 0,4,
                 1,3, 1,5,
                 2,3, 2,6,
                 3,7,
                 4,5, 4,6,
                 5,7,
                 6,7,
                 -1);

    print(&g);

    printf("\nShould be &GY@OQ?OW@?O?:\n");
    igraph_small(&g, 8, IGRAPH_DIRECTED,
                 0,1, 0,2, 0,4,
                 1,3, 1,5,
                 2,3, 2,6,
                 3,7,
                 4,5, 4,6,
                 5,7,
                 6,7,
                 -1);

    print(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
