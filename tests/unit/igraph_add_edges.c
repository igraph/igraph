/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
    igraph_vector_int_t v;

    /* Create graph */
    igraph_small(&g, 0, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 2, 3, 2, 2,
                 -1);

    /* Add edges */
    igraph_vector_int_init_int(&v, 4, 2, 1, 3, 3);
    igraph_add_edges(&g, &v, 0);

    /* Check result */
    igraph_get_edgelist(&g, &v, 0);
    print_vector_int(&v);

    /* Error, vector length */
    igraph_vector_int_resize(&v, 3);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 2;
    CHECK_ERROR(igraph_add_edges(&g, &v, 0), IGRAPH_EINVAL);

    /* Check result */
    igraph_get_edgelist(&g, &v, 0);
    print_vector_int(&v);

    /* Error, vector IDs */
    igraph_vector_int_resize(&v, 4);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 2;
    VECTOR(v)[3] = 4;
    CHECK_ERROR(igraph_add_edges(&g, &v, 0), IGRAPH_EINVVID);

    /* Check result */
    igraph_get_edgelist(&g, &v, 0);
    print_vector_int(&v);

    igraph_vector_int_destroy(&v);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
