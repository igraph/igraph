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

void call_and_print(igraph_t *graph, int size, igraph_vector_t *cut_prob,
                    igraph_integer_t sample_size, igraph_vector_t *parsample) {
    igraph_integer_t estimate;
    IGRAPH_ASSERT(igraph_motifs_randesu_estimate(graph, &estimate, size, cut_prob, sample_size, parsample) == IGRAPH_SUCCESS);
    printf("Estimate: %" IGRAPH_PRId "\n\n", estimate);
}


int main() {
    igraph_t g_0, g_1, g_50_full, g_4_3_1;
    igraph_vector_t cut_prob_0_3;
    igraph_vector_t cut_prob_0_4;
    igraph_vector_t cut_prob_01;
    igraph_vector_t parsample;
    igraph_integer_t estimate;

    igraph_vector_init_real(&cut_prob_0_3, 3, 0.0, 0.0, 0.0);
    igraph_vector_init_real(&cut_prob_0_4, 4, 0.0, 0.0, 0.0, 0.0);
    igraph_vector_init_real(&cut_prob_01, 3, 0.1, 0.1, 0.1);
    igraph_vector_init_seq(&parsample, 0, 40);

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_small(&g_0, 0, 0, -1);
    igraph_small(&g_1, 1, 0, -1);
    igraph_full(&g_50_full, 50, 0, IGRAPH_NO_LOOPS);
    igraph_small(&g_4_3_1, 4, 0, 0,1, 1,2, 2,0, -1);

    printf("No vertices:\n");
    call_and_print(&g_0, /*size*/ 3, &cut_prob_0_3, /*sample_size*/ 1, /*parsample*/ NULL);

    printf("One vertex:\n");
    call_and_print(&g_1, /*size*/ 3, &cut_prob_0_3, /*sample_size*/ 1, /*parsample*/ NULL);

    printf("Full graph of 50 vertices, motif size 3, sample all, (50 choose 3 = 19600):\n");
    call_and_print(&g_50_full, /*size*/ 3, &cut_prob_0_3, /*sample_size*/ 50, /*parsample*/ NULL);

    printf("Full graph of 50 vertices, motif size 3, sample all, cut_prob 0.1 at each level:\n");
    call_and_print(&g_50_full, /*size*/ 3, &cut_prob_01, /*sample_size*/ 50, /*parsample*/ NULL);

    printf("Full graph of 50 vertices, motif size 3, sample 20:\n");
    call_and_print(&g_50_full, /*size*/ 3, &cut_prob_0_3, /*sample_size*/ 20, /*parsample*/ NULL);

    printf("Full graph of 50 vertices, motif size 3, sample first 40:\n");
    call_and_print(&g_50_full, /*size*/ 3, &cut_prob_0_3, /*sample_size*/ 0, &parsample);

    printf("Full graph of 50 vertices, motif size 4, sample 20 (50 choose 4 = 230300:\n");
    call_and_print(&g_50_full, /*size*/ 4, &cut_prob_0_4, /*sample_size*/ 20, /*parsample*/ NULL);

    printf("Triangle and a vertex, motif size 4, sample all:\n");
    call_and_print(&g_4_3_1, /*size*/ 4, &cut_prob_0_4, /*sample_size*/ 4, /*parsample*/ NULL);

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("Cut prob too short.\n");
    IGRAPH_ASSERT(igraph_motifs_randesu_estimate(&g_4_3_1, &estimate, /*size*/ 14, &cut_prob_0_3, /*sample_size*/ 4, /*parsample*/ NULL) == IGRAPH_EINVAL);

    printf("Too many samples.\n");
    IGRAPH_ASSERT(igraph_motifs_randesu_estimate(&g_4_3_1, &estimate, /*size*/ 4, &cut_prob_0_4, /*sample_size*/ 40, /*parsample*/ NULL) == IGRAPH_EINVAL);

    printf("Too many parsamples.\n");
    IGRAPH_ASSERT(igraph_motifs_randesu_estimate(&g_4_3_1, &estimate, /*size*/ 4, &cut_prob_0_4, /*sample_size*/ 4, /*parsample*/ &parsample) == IGRAPH_EINVAL);

    igraph_destroy(&g_0);
    igraph_destroy(&g_1);
    igraph_destroy(&g_50_full);
    igraph_destroy(&g_4_3_1);
    igraph_vector_destroy(&cut_prob_0_3);
    igraph_vector_destroy(&cut_prob_0_4);
    igraph_vector_destroy(&cut_prob_01);
    igraph_vector_destroy(&parsample);

    VERIFY_FINALLY_STACK();
    return 0;
}
