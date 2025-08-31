/*
   igraph library.
   Copyright (C) 2007-2024  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

int main(void) {
    igraph_t graph;
    igraph_vector_int_t membership;
    igraph_real_t modularity;

    /* Initialize the library. */
    igraph_setup();

    igraph_famous(&graph, "Zachary"); /* We use Zachary's karate club network. */

    /* Label propagation is a stochastic method; the result will depend on the random seed. */
    igraph_rng_seed(igraph_rng_default(), 123);

    /* All igraph functions that returns their result in an igraph_vector_t must be given
       an already initialized vector. */
    igraph_vector_int_init(&membership, 0);
    igraph_community_label_propagation(
        &graph, &membership, /* mode = */ IGRAPH_ALL,
        /* weights= */ NULL, /* initial= */ NULL, /* fixed= */ NULL,
        IGRAPH_LPA_DOMINANCE);

    /* Also calculate the modularity of the partition */
    igraph_modularity(
        &graph, &membership, /* weights= */ NULL, /* resolution = */ 1,
        /* directed= */ false, &modularity);

    printf("%" IGRAPH_PRId " communities found; modularity score is %g.\n",
           igraph_vector_int_max(&membership) + 1, modularity);

    printf("Communities membership: ");
    igraph_vector_int_print(&membership);

    /* Destroy data structures at the end. */
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&graph);

    return 0;
}
