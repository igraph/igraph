/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2022  The igraph development team

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

int main(void) {
    igraph_t g;
    igraph_vector_t degree;
    igraph_plfit_result_t model;

    /* Seed random number generator to ensure reproducibility. */
    igraph_rng_seed(igraph_rng_default(), 42);

    /* Generate a BA network; degree distribution is supposed to be a power-law
     * if the graph is large enough */
    igraph_barabasi_game(
        &g, 10000, /*power=*/ 1, /*m=*/ 2,
        /* outseq= */ 0, /* outpref= */ 0, /*A=*/ 1,
        IGRAPH_UNDIRECTED, IGRAPH_BARABASI_BAG,
        /*start_from=*/ 0
    );

    /* Get the vertex degrees. We use igraph_strength() because it stores its
     * result in an igraph_vector_t */
    igraph_vector_init(&degree, 0);
    igraph_strength(&g, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS, 0);

    /* Fit a power-law to the degrees */
    igraph_power_law_fit(
        &degree, &model, /* xmin = */ -1,
        /* force_continuous = */ 0
    );

    /* If you also need a p-value: */
    /* igraph_plfit_result_calculate_p_value(&model, &p, 0.001); */

    printf("alpha = %.5f\n", model.alpha);
    printf("xmin = %.5f\n", model.xmin);
    printf("log-likelihood = %.5f\n", model.L);

    igraph_vector_destroy(&degree);
    igraph_destroy(&g);

    return 0;
}
