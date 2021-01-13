/*
   IGraph library.
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
#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_strvector_t names, weights;

    /*
    Expected output:

    ```
    # vertex1name
    vertex2name [optionalWeight]
    vertex3name [optionalWeight]
    ```

    where vertex 2 and 3 are connected with vectex 1 by an edge.

    Isolated vertices are optionally shown.
    */

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_small(&g, 7, IGRAPH_UNDIRECTED, 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 3, 4, -1);


    printf("Test without isolates:\n");

    igraph_write_graph_lgl(&g, stdout, /*names*/ NULL, /*weights*/ NULL, /*isolates*/ 0);


    printf("\nTest with isolates:\n");

    igraph_write_graph_lgl(&g, stdout, /*names*/ NULL, /*weights*/ NULL, /*isolates*/ 1);


    printf("\nTest vertex and edge labels:\n");

    igraph_strvector_init(&names, 7);
    igraph_strvector_set(&names, 0, "zero");
    igraph_strvector_set(&names, 1, "one");
    igraph_strvector_set(&names, 2, "two");
    igraph_strvector_set(&names, 3, "three");
    igraph_strvector_set(&names, 4, "four");
    igraph_strvector_set(&names, 5, "five");
    igraph_strvector_set(&names, 6, "six");
    SETVASV(&g, "names", &names);

    igraph_strvector_init(&weights, 6);
    igraph_strvector_set(&weights, 0, "8");
    igraph_strvector_set(&weights, 1, "9");
    igraph_strvector_set(&weights, 2, "10");
    igraph_strvector_set(&weights, 3, "11");
    igraph_strvector_set(&weights, 4, "12");
    igraph_strvector_set(&weights, 5, "13");

    SETEASV(&g, "weights", &weights);

    igraph_write_graph_lgl(&g, stdout, "names", "weights", /*isolates*/ 0);

    igraph_strvector_destroy(&names);
    igraph_strvector_destroy(&weights);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
