/*
   igraph library.
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

void compute_and_print(const igraph_t *g, igraph_connectedness_t mode) {
    igraph_vector_int_t membership, csize;
    igraph_int_t no;

    igraph_vector_int_init(&membership, 1);
    igraph_vector_int_init(&csize, 1);

    igraph_connected_components(g, &membership, &csize, &no, mode);

    /* Canonicalize membership representation */
    igraph_reindex_membership(&membership, NULL, NULL);

    printf("Mode: %s\n", mode == IGRAPH_STRONG ? "strong" : "weak");
    print_vector_int(&membership);
    print_vector_int(&csize);
    printf("No. of components: %" IGRAPH_PRId "\n", no);

    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&membership);
}

int main(void) {
    igraph_t g;

    /* Component calculations are run twice to exercise the cache */

    printf("\nNull graph (not connected)\n");
    igraph_empty(&g, 0, IGRAPH_DIRECTED);
    compute_and_print(&g, IGRAPH_STRONG);
    compute_and_print(&g, IGRAPH_STRONG);
    igraph_invalidate_cache(&g);
    compute_and_print(&g, IGRAPH_WEAK);
    compute_and_print(&g, IGRAPH_WEAK);
    igraph_destroy(&g);

    printf("\nSingleton graph (connected)\n");
    igraph_empty(&g, 1, IGRAPH_DIRECTED);
    compute_and_print(&g, IGRAPH_STRONG);
    compute_and_print(&g, IGRAPH_STRONG);
    igraph_invalidate_cache(&g);
    compute_and_print(&g, IGRAPH_WEAK);
    compute_and_print(&g, IGRAPH_WEAK);
    igraph_destroy(&g);

    printf("\nKautz graph (connected)\n");
    igraph_kautz(&g, 2, 2);
    compute_and_print(&g, IGRAPH_STRONG);
    compute_and_print(&g, IGRAPH_STRONG);
    igraph_invalidate_cache(&g);
    compute_and_print(&g, IGRAPH_WEAK);
    compute_and_print(&g, IGRAPH_WEAK);
    igraph_destroy(&g);

    printf("\nDirected 2-path\n");
    igraph_small(&g, 2, IGRAPH_DIRECTED, 0,1, -1);
    compute_and_print(&g, IGRAPH_STRONG);
    compute_and_print(&g, IGRAPH_STRONG);
    igraph_invalidate_cache(&g);
    compute_and_print(&g, IGRAPH_WEAK);
    compute_and_print(&g, IGRAPH_WEAK);
    igraph_destroy(&g);

    printf("\nTwo disjoint 3-cycles\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, 3,4, 4,5, 5,3, -1);
    compute_and_print(&g, IGRAPH_STRONG);
    compute_and_print(&g, IGRAPH_STRONG);
    igraph_invalidate_cache(&g);
    compute_and_print(&g, IGRAPH_WEAK);
    compute_and_print(&g, IGRAPH_WEAK);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
