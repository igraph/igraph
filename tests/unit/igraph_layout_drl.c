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

#include <igraph.h>
#include "test_utilities.inc"

void set_options_fast(igraph_layout_drl_options_t *options) {
        options->edge_cut = 4.0/5.0;

        options->init_iterations   = 10;
        options->init_temperature  = 2000;
        options->init_attraction   = 10;
        options->init_damping_mult = 1.0;

        options->liquid_iterations   = 10;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 10;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = 10;
        options->expansion_temperature  = 2000;
        options->expansion_attraction   = 2;
        options->expansion_damping_mult = 1.0;

        options->cooldown_iterations   = 10;
        options->cooldown_temperature  = 2000;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = 10;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = 10;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = 1;
}

void check_and_destroy(igraph_matrix_t *result, igraph_real_t half_size) {
    igraph_real_t min, max;
    igraph_matrix_minmax(result, &min, &max);
    IGRAPH_ASSERT(min >= -half_size);
    IGRAPH_ASSERT(max <= half_size);
    igraph_matrix_destroy(result);
}

int main() {
    igraph_t g;
    igraph_matrix_t result;
    igraph_layout_drl_options_t options;
    int i;
    igraph_real_t *damping_muls[6] = {&options.init_damping_mult, &options.liquid_damping_mult, &options.expansion_damping_mult, &options.cooldown_damping_mult, &options.crunch_damping_mult, &options.simmer_damping_mult};

    igraph_rng_seed(igraph_rng_default(), 42);

    set_options_fast(&options);

    printf("The Zachary karate club.\n");
    igraph_famous(&g, "zachary");
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 50);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Negative damping.\n");
    igraph_matrix_init(&result, 0, 0);
    for (i = 0; i < 6; i++) {
        *damping_muls[i] *= -1.0;
        IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                    /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_EINVAL);
        *damping_muls[i] *= -1.0;
    }
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
