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

int main(void) {
    igraph_t g, gzero, h;
    igraph_vector_t quant, quantzero;
    igraph_vector_int_t strat, stratzero;
    igraph_integer_t nvert;

    /* nonempty graph */
    igraph_small(&g, /*nvert=*/ 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_empty(&h, 0, 0);         /* empty graph */
    igraph_vector_init(&quant, 1);  /* quantities vector */
    igraph_vector_int_init(&strat, 2);  /* strategies vector */
    igraph_small(&gzero, /*nvert=*/ 0, IGRAPH_UNDIRECTED,
                 0, 3, 0, 4, 1, 2, 1, 4, 1, 5, 2, 3, 2, 4, 3, 4, -1);
    nvert = igraph_vcount(&gzero);
    igraph_vector_int_init_int(&stratzero, nvert, 1, 0, 1, 2, 0, 3);
    igraph_vector_init(&quantzero, nvert);  /* vector of zeros */

    /* test parameters */
    /*graph--vert--islocal--quantities--strategies--mode--retval*/
    /* null pointer for graph */
    CHECK_ERROR(igraph_roulette_wheel_imitation(NULL, 0, 1, NULL, NULL, IGRAPH_ALL), IGRAPH_EINVAL);
    /* null pointer for quantities vector */
    CHECK_ERROR(igraph_roulette_wheel_imitation(&g, 0, 1, NULL, NULL, IGRAPH_ALL), IGRAPH_EINVAL);
    /* null pointer for strategies vector */
    CHECK_ERROR(igraph_roulette_wheel_imitation(&g, 0, 1, &quant, NULL, IGRAPH_ALL), IGRAPH_EINVAL);
    /* empty graph */
    CHECK_ERROR(igraph_roulette_wheel_imitation(&h, 0, 1, &quant, &strat, IGRAPH_ALL), IGRAPH_EINVAL);
    /* length of quantities vector different from number of vertices */
    CHECK_ERROR(igraph_roulette_wheel_imitation(&g, 0, 1, &quant, &strat, IGRAPH_ALL), IGRAPH_EINVAL);
    /* length of strategies vector different from number of vertices */
    CHECK_ERROR(igraph_roulette_wheel_imitation(&g, 0, 1, &quant, &strat, IGRAPH_ALL), IGRAPH_EINVAL);
    /* quantities vector contains all zeros */
    CHECK_ERROR(igraph_roulette_wheel_imitation(&gzero, 4, 1, &quantzero, &stratzero, IGRAPH_ALL), IGRAPH_EINVAL);

    igraph_destroy(&g);
    igraph_destroy(&gzero);
    igraph_destroy(&h);
    igraph_vector_destroy(&quant);
    igraph_vector_destroy(&quantzero);
    igraph_vector_int_destroy(&strat);
    igraph_vector_int_destroy(&stratzero);

    VERIFY_FINALLY_STACK();

    return 0;
}
