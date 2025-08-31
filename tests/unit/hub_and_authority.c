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


void print_hub_and_authority(igraph_t *g, igraph_vector_t *weights, igraph_bool_t use_options) {
    igraph_arpack_options_t options;
    igraph_vector_t hub_vector, authority_vector;
    igraph_real_t value;

    igraph_arpack_options_init(&options);
    igraph_vector_init(&hub_vector, 0);
    igraph_vector_init(&authority_vector, 0);

    printf("--------------------------------------------------\n");

    igraph_hub_and_authority_scores(g, &hub_vector, &authority_vector, &value,
            weights, use_options ? &options : NULL);

    vector_chop(&hub_vector, 10e-10);
    vector_chop(&authority_vector, 10e-10);
    printf("hub:\n");
    print_vector(&hub_vector);
    printf("authority:\n");
    print_vector(&authority_vector);
    printf("value:\n");
    print_real(stdout, value, "%g");
    printf("\n");
    printf("--------------------------------------------------\n\n\n");
    igraph_vector_destroy(&hub_vector);
    igraph_vector_destroy(&authority_vector);
}

int main(void) {
    igraph_t g;
    igraph_vector_t weights;
    igraph_arpack_options_t options;
    igraph_real_t value;

    igraph_arpack_options_init(&options);

    printf("Null graph:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
    print_hub_and_authority(&g, NULL, true);
    igraph_destroy(&g);

    printf("Singleton graph with loop:\n");
    igraph_small(&g, 1, IGRAPH_DIRECTED, 0,0, -1);
    print_hub_and_authority(&g, NULL, true);
    igraph_destroy(&g);

    printf("Singleton graph with three loops:\n");
    igraph_small(&g, 1, IGRAPH_DIRECTED, 0,0, 0,0, 0,0, -1);
    print_hub_and_authority(&g, NULL, true);
    igraph_destroy(&g);

    printf("Three vertices, no links:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED, -1);
    print_hub_and_authority(&g, NULL, true);
    igraph_destroy(&g);

    printf("Two hubs and one authority:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED,
        0,2, 1,2, -1);
    igraph_vector_init_int(&weights, 2,
        1, 1);
    print_hub_and_authority(&g, &weights, true);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    /* https://nlp.stanford.edu/IR-book/html/htmledition/hubs-and-authorities-1.html  */
    /* with different normalization */
    printf("Stanford example:\n");
    igraph_small(&g, 7, IGRAPH_DIRECTED,
        0,2, 1,1, 1,2, 2,0, 2,2, 2,3, 3,3, 3,4, 4,6, 5,5,
        5,6, 6,3, 6,4, 6,6, -1);
    igraph_vector_init_int(&weights, 14,
        1, 1, 1, 1, 1, 2, 1, 1, 1, 1,
        1, 2, 1, 1);
    print_hub_and_authority(&g, &weights, false);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    printf("Same example with scaling:\n");
    igraph_small(&g, 7, IGRAPH_DIRECTED,
        0,2, 1,1, 1,2, 2,0, 2,2, 2,3, 3,3, 3,4, 4,6, 5,5,
        5,6, 6,3, 6,4, 6,6, -1);
    igraph_vector_init_int(&weights, 14,
        1, 1, 1, 1, 1, 2, 1, 1, 1, 1,
        1, 2, 1, 1);
    print_hub_and_authority(&g, &weights, false);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    /* Verify that self-loops are counted twice in undirected graphs. */
    printf("Undirected graph with self-loops and multi-edges:\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 2, 2, 3, 2, 3,
                 -1);
    print_hub_and_authority(&g, NULL, false);
    igraph_destroy(&g);

    printf("Degenerate example:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED,
        0,1, 1,0, 1,2, 2,1, 2,3, 3,0, -1);
    igraph_hub_and_authority_scores(&g, NULL, NULL, &value,
                                    NULL, &options);
    printf("--------------------------------------------------\n");
    printf("value:\n");
    print_real(stdout, value, "%g");
    printf("\n");
    printf("--------------------------------------------------\n\n\n");
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Checking invalid weight vector.\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED,
        0,2, 1,2, -1);
    igraph_vector_init_int(&weights, 3,
        1, 1, 1);
    IGRAPH_ASSERT(igraph_hub_and_authority_scores(&g, NULL, NULL, NULL,
                                    &weights, &options) == IGRAPH_EINVAL);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    VERIFY_FINALLY_STACK();
    return 0;
}
