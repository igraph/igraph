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

void check_and_destroy(igraph_matrix_t *result, igraph_real_t half_size) {
    igraph_real_t min, max;
    igraph_matrix_minmax(result, &min, &max);
    IGRAPH_ASSERT(min >= -half_size);
    IGRAPH_ASSERT(max <= half_size);
    igraph_matrix_destroy(result);
}

int main(void) {
    igraph_t g;
    igraph_matrix_t result;
    igraph_layout_drl_options_t options;
    int i;
    igraph_real_t *damping_muls[6] = {&options.init_damping_mult, &options.liquid_damping_mult, &options.expansion_damping_mult, &options.cooldown_damping_mult, &options.crunch_damping_mult, &options.simmer_damping_mult};

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_REFINE);

    printf("The Zachary karate club.\n");
    igraph_famous(&g, "zachary");
    igraph_matrix_init(&result, 0, 0);
    igraph_layout_drl_3d(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL);
    check_and_destroy(&result, 50);

    VERIFY_FINALLY_STACK();

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_FINAL);
    printf("Negative damping.\n");
    igraph_matrix_init(&result, 0, 0);
    for (i = 0; i < 6; i++) {
        *damping_muls[i] = -1.0;
        CHECK_ERROR(igraph_layout_drl_3d(&g, &result, /*use_seed*/ 0, &options,
                    /*weights*/ NULL), IGRAPH_EINVAL);
        *damping_muls[i] = 1.0;
    }
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
