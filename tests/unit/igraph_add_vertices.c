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

    igraph_t g1;
    igraph_vector_int_t v1;

    /* Create a graph */
    igraph_vector_int_init(&v1, 8);
    VECTOR(v1)[0] = 0;
    VECTOR(v1)[1] = 1;
    VECTOR(v1)[2] = 1;
    VECTOR(v1)[3] = 2;
    VECTOR(v1)[4] = 2;
    VECTOR(v1)[5] = 3;
    VECTOR(v1)[6] = 2;
    VECTOR(v1)[7] = 2;
    igraph_create(&g1, &v1, 0, 0);
    igraph_vector_int_destroy(&v1);

    /* Add more vertices */
    igraph_add_vertices(&g1, 10, 0);
    IGRAPH_ASSERT(igraph_vcount(&g1) == 14);

    /* Add more vertices */
    igraph_add_vertices(&g1, 0, 0);
    IGRAPH_ASSERT(igraph_vcount(&g1) == 14);

    /* Error */
    CHECK_ERROR(igraph_add_vertices(&g1, -1, 0), IGRAPH_EINVAL);

    igraph_destroy(&g1);

    VERIFY_FINALLY_STACK();

    return 0;
}
