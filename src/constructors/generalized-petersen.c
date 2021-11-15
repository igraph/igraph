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

#include "igraph_constructors.h"

#include "igraph_interface.h"


int igraph_generalized_petersen(igraph_t *igraph, igraph_integer_t n, igraph_integer_t k) {
    /* This is a generalized Petersen graph constructor*/
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (n < 3) {
        IGRAPH_ERROR("n must be greater than or equal to 3", IGRAPH_EINVAL);
    }

    if (k < 1 || k > (n-1)/2) {
        IGRAPH_ERROR("k must be positive and less than (n-1)/2", IGRAPH_EINVAL);
    }




    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    /* Return no error */
    return IGRAPH_SUCCESS;
}
