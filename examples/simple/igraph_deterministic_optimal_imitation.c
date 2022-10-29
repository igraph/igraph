/* -*- mode: C -*-  */
/*
  Test suite for deterministic optimal imitation.
  Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>

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
#include <stdio.h>

/* test parameters structure */
typedef struct {
    igraph_t *graph;
    igraph_integer_t vertex;
    igraph_optimal_t optimality;
    igraph_vector_t *quantities;
    igraph_vector_int_t *strategies;
    igraph_neimode_t mode;
    igraph_error_t retval;
} strategy_test_t;

/* Updating the strategy of an isolated vertex. In this case, the strategies
 * vector should not change at all.
 */
igraph_error_t isolated_vertex_test(void) {
    igraph_t g;
    igraph_vector_t quant;
    igraph_vector_int_t strat, v;
    igraph_integer_t i;
    igraph_error_t ret;

    /* graph with one isolated vertex */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_add_vertices(&g, 1, 0);  /* new vertex 3 is isolated */
    /* quantities vector: all vertices have the same fitness */
    igraph_vector_init_real(&quant, 4, 0.25, 0.25, 0.25, 0.25);
    /* strategies vector: 0 means aggressive strategy; 1 means passive */
    igraph_vector_int_init_int(&strat, 4, 1, 0, 1, 0);
    /* make a copy of the original strategies vector for comparison later on */
    igraph_vector_int_init_copy(&v, &strat);
    /* Now update strategy of vertex 3. Since this vertex is isolated, no */
    /* strategy update would take place. The resulting strategies vector */
    /* would be the same as it was originally. */
    ret = igraph_deterministic_optimal_imitation(/*graph*/ &g,
            /*vertex*/ 3,
            /*optimality*/ IGRAPH_MAXIMUM,
            /*quantities*/ &quant,
            /*strategies*/ &strat,
            /*mode*/ IGRAPH_ALL);
    if (ret) {
        printf("Isolated vertex test failed.\n");
        return IGRAPH_FAILURE;
    }
    for (i = 0; i < igraph_vector_int_size(&strat); i++) {
        if (VECTOR(strat)[i] != VECTOR(v)[i]) {
            printf("Isolated vertex test failed.\n");
            return IGRAPH_FAILURE;
        }
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);
    igraph_vector_int_destroy(&v);

    return IGRAPH_SUCCESS;
}

/* A game on the Petersen graph. This graph has 10 vertices and 15 edges.
 * The Petersen graph is initialized with a default quantities vector and a
 * default strategies vector. For each vertex v in the graph, we update the
 * strategy of v via deterministic optimal imitation. The resulting updated
 * strategies vector is compared with the known result vector. A mismatch would
 * raise an error code. If the updated strategies vector matches the known
 * result vector, we reset the strategies vector to its default state and
 * repeat the game with another vertex.
 */
igraph_error_t petersen_game_test(void) {
    igraph_t g;
    igraph_vector_t quant;
    igraph_vector_int_t known_max_v, known_min_v, strat, stratcopy;
    igraph_integer_t i, nvert;

    /* the Petersen graph */
    igraph_small(&g, /*n=*/ 0, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 5, 1, 2, 1, 6, 2, 3, 2, 7, 3, 4, 3, 8, 4, 9,
                 5, 7, 5, 8, 6, 8, 6, 9, 7, 9, -1);
    nvert = igraph_vcount(&g);
    /* Strategies vector, one strategy for each vertex. Thus vec[i] is the */
    /* strategy of vertex i. The strategy space is: {0, 1, 2, 3}. */
    igraph_vector_int_init_int(&strat, nvert, 1, 1, 2, 2, 0, 0, 0, 1, 2, 3);
    /* Quantities vector, one quantity per vertex. Thus vec[i] is the */
    /* quantity for vertex i. */
    igraph_vector_init_real(&quant, nvert,
                            0.3, 1.1, 0.5, 1.0, 0.9,
                            0.8, 0.4, 0.1, 0.7, 0.7);
    /* Known strategies that would be adopted. Thus vec[i] means that in */
    /* game i where we revise the strategy of vertex i, the strategy */
    /* vec[i] would be adopted by i. */
    /*maximum deterministic imitation*/
    igraph_vector_int_init_int(&known_max_v, nvert, 1, 1, 1, 2, 2, 0, 1, 0, 2, 0);
    /*minimum deterministic imitation*/
    igraph_vector_int_init_int(&known_min_v, nvert, 1, 1, 1, 2, 1, 1, 0, 1, 0, 1);
    /* play game and compare resulting updated strategies */
    for (i = 0; i < nvert; i++) {
        /* maximum deterministic imitation */
        igraph_vector_int_init_copy(&stratcopy, &strat);
        igraph_deterministic_optimal_imitation(/*graph*/ &g,
                /*vertex*/ i,
                /*optimality*/ IGRAPH_MAXIMUM,
                /*quantities*/ &quant,
                /*strategies*/ &stratcopy,
                /*neighbours*/ IGRAPH_ALL);
        if (VECTOR(stratcopy)[i] != VECTOR(known_max_v)[i]) {
            printf("Maximum deterministic imitation failed for vertex %" IGRAPH_PRId ".\n", i);
            return IGRAPH_FAILURE;
        }
        igraph_vector_int_destroy(&stratcopy);
        /* minimum deterministic imitation */
        igraph_vector_int_init_copy(&stratcopy, &strat);
        igraph_deterministic_optimal_imitation(/*graph*/ &g,
                /*vertex*/ i,
                /*optimality*/ IGRAPH_MINIMUM,
                /*quantities*/ &quant,
                /*strategies*/ &stratcopy,
                /*neighbours*/ IGRAPH_ALL);
        if (VECTOR(stratcopy)[i] != VECTOR(known_min_v)[i]) {
            printf("Minimum deterministic imitation failed for vertex %" IGRAPH_PRId ".\n", i);
            return IGRAPH_FAILURE;
        }
        igraph_vector_int_destroy(&stratcopy);
    }
    /* clean up */
    igraph_destroy(&g);
    igraph_vector_int_destroy(&known_max_v);
    igraph_vector_int_destroy(&known_min_v);
    igraph_vector_destroy(&quant);
    igraph_vector_int_destroy(&strat);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_error_t ret;

    igraph_rng_seed(igraph_rng_default(), 648);

    ret = isolated_vertex_test();
    if (ret) {
        return ret;
    }
    ret = petersen_game_test();
    if (ret) {
        return ret;
    }

    return IGRAPH_SUCCESS;
}
