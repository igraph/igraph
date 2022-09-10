/* IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

int main() {
    igraph_vector_int_t v, seq;
    igraph_t g;

    igraph_vector_int_init(&seq, 3);
    igraph_vector_int_init(&v, 0);
    VECTOR(seq)[0] = 2;
    VECTOR(seq)[1] = 0;
    VECTOR(seq)[2] = 2;

    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 2,2, -1);
    /* Invalid mode */
    CHECK_ERROR(igraph_degree(&g, &v, igraph_vss_vector(&seq), (igraph_neimode_t)0,
                        IGRAPH_LOOPS), IGRAPH_EINVMODE);

    VECTOR(seq)[0] = 4;
    /* Vertex does not exist */
    CHECK_ERROR(igraph_degree(&g, &v, igraph_vss_vector(&seq), IGRAPH_ALL, IGRAPH_LOOPS), IGRAPH_EINVVID);

    igraph_destroy(&g);
    igraph_vector_int_destroy(&v);
    igraph_vector_int_destroy(&seq);

    VERIFY_FINALLY_STACK();

    return 0;
}
