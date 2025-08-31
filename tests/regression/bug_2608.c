
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

#include "../unit/test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t membership;

    /* label propagation is a stochastic method */
    igraph_rng_seed(igraph_rng_default(), 765);

    igraph_small(&g, 0, /* directed = */ 1,
        0, 1, 0, 2, 1, 0, 1, 0, 1, 1, 2, 0, 2, 1, 2, 1, 2, 1, -1);
    igraph_vector_int_init(&membership, 0);

    igraph_community_label_propagation(&g, &membership, IGRAPH_OUT, NULL, NULL, NULL, IGRAPH_LPA_FAST);

    print_vector_int(&membership);
    igraph_vector_int_destroy(&membership);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
