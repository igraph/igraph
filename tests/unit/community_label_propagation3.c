/*
   igraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t membership, initial;
    igraph_vector_bool_t fixed;
    igraph_int_t i;

    /* label propagation is a stochastic method */
    igraph_rng_seed(igraph_rng_default(), 765);

    /* Zachary Karate club -- does not matter, we are simply interested in
     * whether the system handles it correctly if a fixed node is unlabeled */
    igraph_famous(&g, "zachary");


    igraph_vector_int_init(&initial, igraph_vcount(&g));
    igraph_vector_int_fill(&initial, -1);
    igraph_vector_bool_init(&fixed, igraph_vcount(&g));
    igraph_vector_bool_fill(&fixed, 0);
    VECTOR(fixed)[7] = 1;
    VECTOR(fixed)[13] = 1;

    igraph_vector_int_init(&membership, 0);

    igraph_community_label_propagation(&g, &membership, IGRAPH_OUT, NULL, &initial, &fixed, IGRAPH_LPA_DOMINANCE);

    for (i = 0; i < igraph_vcount(&g); i++) {
        /* Check that the "fixed" vector has not been changed */
        if (i == 7 || i == 13) {
            IGRAPH_ASSERT(VECTOR(fixed)[i]);
        } else {
            IGRAPH_ASSERT(!VECTOR(fixed)[i]);
        }

        /* Check that no vertex remained unlabeled */
        IGRAPH_ASSERT(VECTOR(membership)[i] >= 0);
    }

    igraph_vector_bool_destroy(&fixed);
    igraph_vector_int_destroy(&initial);
    igraph_vector_int_destroy(&membership);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
