/* igraph library.
   Copyright (C) 2010-2024  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t graph;
    igraph_vector_int_t membership;
    igraph_real_t modularity;
    igraph_error_handler_t *handler;
    igraph_error_t errcode;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&graph, 9, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 0, 1, 0, 2, 0, 5, 0, 6, 1, 3, 3, 5, 4, 7, 7, 8, 2, 4, 1, 2, 6, 8,
                 -1);

    igraph_vector_int_init(&membership, 0);

    /* Temporarily ignore errors so we can check if the functionality is available. */
    handler = igraph_set_error_handler(&igraph_error_handler_printignore);

    errcode = igraph_community_optimal_modularity(&graph, NULL, 1, &modularity, &membership);
    if (errcode == IGRAPH_UNIMPLEMENTED) {
        /* igraph was not build with GLPK; functionality unavailable, so we quit */
        return 77;
    }

    /* Restore the previous error handler. */
    igraph_set_error_handler(handler);

    printf("Highest possible modularity: %g\n", modularity);
    printf("Membership: ");
    igraph_vector_int_print(&membership);

    igraph_destroy(&graph);

    igraph_vector_int_destroy(&membership);

    return 0;
}
