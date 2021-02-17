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

void set_options(igraph_layout_drl_options_t *options) {
        options->init_iterations   = -0;
        options->init_temperature  = -2000;
        options->init_attraction   = -10;
        options->init_damping_mult = 1.0;

        options->liquid_iterations   = -200;
        options->liquid_temperature  = -2000;
        options->liquid_attraction   = -10;
        options->liquid_damping_mult = 1.0;

        options->expansion_iterations   = -200;
        options->expansion_temperature  = -2000;
        options->expansion_attraction   = -2;
        options->expansion_damping_mult = 1.0;

        options->cooldown_iterations   = -200;
        options->cooldown_temperature  = -2000;
        options->cooldown_attraction   = -1;
        options->cooldown_damping_mult = .1;

        options->crunch_iterations   = -50;
        options->crunch_temperature  = -250;
        options->crunch_attraction   = -1;
        options->crunch_damping_mult = 0.25;

        options->simmer_iterations   = -100;
        options->simmer_temperature  = -250;
        options->simmer_attraction   = -.5;
        options->simmer_damping_mult = 0;
}

void set_options_2(igraph_layout_drl_options_t *options) {
        options->init_iterations   = 0;
        options->init_temperature  = 2000;
        options->init_attraction   = 10;
        options->init_damping_mult = -1.0;

        options->liquid_iterations   = 200;
        options->liquid_temperature  = 2000;
        options->liquid_attraction   = 10;
        options->liquid_damping_mult = -1.0;

        options->expansion_iterations   = 200;
        options->expansion_temperature  = 2000;
        options->expansion_attraction   = 2;
        options->expansion_damping_mult = -1.0;

        options->cooldown_iterations   = 200;
        options->cooldown_temperature  = 2000;
        options->cooldown_attraction   = 1;
        options->cooldown_damping_mult = -.1;

        options->crunch_iterations   = 50;
        options->crunch_temperature  = 250;
        options->crunch_attraction   = 1;
        options->crunch_damping_mult = -0.25;

        options->simmer_iterations   = 100;
        options->simmer_temperature  = 250;
        options->simmer_attraction   = .5;
        options->simmer_damping_mult = -0;
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
    igraph_vector_t weights;
    igraph_layout_drl_options_t options;
    igraph_real_t seed[20] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0};
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_DEFAULT);

    printf("Empty graph.\n");
    igraph_small(&g, 0, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_SUCCESS);
    print_matrix(&result);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    printf("1 vertex.\n");
    igraph_small(&g, 1, 0, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 100);
    igraph_destroy(&g);

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_REFINE);

    printf("A disconnected graph of 10 vertices with loops.\n");
    igraph_small(&g, 10, 0, 0,1, 1,2, 2,0, 5,6, 6,7, 7,6, 7,7, 8,8, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 500);
    igraph_destroy(&g);

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_COARSEN);

    printf("A disconnected graph of 10 vertices with loops and weights and seed.\n");
    igraph_vector_init(&weights, 8);
    igraph_vector_fill(&weights, 100);
    igraph_small(&g, 10, 0, 0,1, 1,2, 2,0, 5,6, 6,7, 7,6, 7,7, 8,8, -1);
    matrix_init_real_row_major(&result, 10, 2, seed);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 1, &options,
                  &weights, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 500);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_FINAL);

    printf("A disconnected graph of 10 vertices with loops and weights and seed with different options.\n");
    igraph_vector_init(&weights, 8);
    igraph_vector_fill(&weights, 100);
    igraph_small(&g, 10, 0, 0,1, 1,2, 2,0, 5,6, 6,7, 7,6, 7,7, 8,8, -1);
    matrix_init_real_row_major(&result, 10, 2, seed);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 1, &options,
                  &weights, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 500);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    igraph_layout_drl_options_init(&options, IGRAPH_LAYOUT_DRL_COARSEST);

    printf("Negative weights.\n");
    igraph_vector_init(&weights, 8);
    igraph_vector_fill(&weights, -100);
    igraph_small(&g, 10, 0, 0,1, 1,2, 2,0, 5,6, 6,7, 7,6, 7,7, 8,8, -1);
    igraph_matrix_init(&result, 0, 0);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  &weights, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 2000);
    igraph_destroy(&g);
    igraph_vector_destroy(&weights);

    printf("Negative options, positive damping.\n");
    igraph_small(&g, 10, 0, 0,1, 1,2, 2,0, 5,6, 6,7, 7,6, 7,7, 8,8, -1);
    igraph_matrix_init(&result, 0, 0);
    set_options(&options);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_SUCCESS);
    check_and_destroy(&result, 50);
    igraph_destroy(&g);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Negative damping.\n");
    igraph_small(&g, 10, 0, 0,1, 1,2, 2,0, 5,6, 6,7, 7,6, 7,7, 8,8, -1);
    igraph_matrix_init(&result, 0, 0);
    set_options_2(&options);
    IGRAPH_ASSERT(igraph_layout_drl(&g, &result, /*use_seed*/ 0, &options,
                  /*weights*/ NULL, /*fixed*/ 0) == IGRAPH_EINVAL);
    igraph_matrix_destroy(&result);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
