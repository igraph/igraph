/*
   IGraph library.
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

int main(void) {
    igraph_t g, h;
    igraph_vector_t quant;
    igraph_vector_int_t strat;

    /* nonempty graph */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_empty(&h, 0, 0);         /* empty graph */
    igraph_vector_init(&quant, 1);  /* quantities vector */
    igraph_vector_int_init(&strat, 2);  /* strategies vector */

    /* test parameters */
    /*--graph--vertex--optimality--quantities--strategies--mode--retval--*/
    /* null pointer for graph */
    CHECK_ERROR(igraph_deterministic_optimal_imitation( NULL, 0, IGRAPH_MINIMUM, NULL, NULL, IGRAPH_ALL), IGRAPH_EINVAL );
    /* null pointer for quantities vector */
    CHECK_ERROR(igraph_deterministic_optimal_imitation( &g, 0, IGRAPH_MINIMUM, NULL, NULL, IGRAPH_ALL), IGRAPH_EINVAL );
    /* null pointer for strategies vector */
    CHECK_ERROR(igraph_deterministic_optimal_imitation( &g, 0, IGRAPH_MINIMUM, &quant, NULL, IGRAPH_ALL), IGRAPH_EINVAL );
    /* empty graph */
    CHECK_ERROR(igraph_deterministic_optimal_imitation(&h, 0, IGRAPH_MINIMUM, &quant, &strat, IGRAPH_ALL), IGRAPH_EINVAL );
    /* length of quantities vector different from number of vertices */
    CHECK_ERROR(igraph_deterministic_optimal_imitation(&g, 0, IGRAPH_MINIMUM, &quant, &strat, IGRAPH_ALL), IGRAPH_EINVAL );
    /* length of strategies vector different from number of vertices */
    CHECK_ERROR(igraph_deterministic_optimal_imitation(&g, 0, IGRAPH_MINIMUM, &quant, &strat, IGRAPH_ALL), IGRAPH_EINVAL );

    /* clean up */
    igraph_destroy(&g);
    igraph_destroy(&h);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);

    VERIFY_FINALLY_STACK();

    return 0;
}
