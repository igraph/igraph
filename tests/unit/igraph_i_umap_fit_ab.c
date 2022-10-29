/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <stdlib.h>
#include <math.h>

#include "layout/layout_internal.h"

#include "test_utilities.h"

int main(void) {
    size_t i;
    igraph_real_t a, b;

    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */
    igraph_real_t min_dists[8] = {0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0};

    RNG_BEGIN();

    /* test with various typical min_dist values. Originally there is a scaling sigma
     * factor, but it's 1.0 in all default cases so we fix it for now */
    for (i = 0; i < sizeof(min_dists) / sizeof(min_dists[0]); i++) {
        igraph_i_umap_fit_ab(min_dists[i], &a, &b);
        printf("%g, %.1g, %.1g\n", min_dists[i], a, b);
    }

    RNG_END();

    VERIFY_FINALLY_STACK();

    return 0;
}
