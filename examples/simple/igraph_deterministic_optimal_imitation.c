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
  igraph_vector_t *strategies;
  igraph_neimode_t mode;
} strategy_test_t;

/* Error tests. That is, we expect error codes to be returned from such tests.
 */
int error_tests() {
  igraph_t g, h;
  igraph_vector_t quant, strat,  v;
  int i, n, ret;
  strategy_test_t *test;

  igraph_vector_init(&v, 6);
  VECTOR(v)[0] = 0; VECTOR(v)[1] = 1;
  VECTOR(v)[2] = 1; VECTOR(v)[3] = 2;
  VECTOR(v)[4] = 2; VECTOR(v)[5] = 0;
  igraph_create(&g, &v, 0, 0);    /* nonempty graph */
  igraph_vector_destroy(&v);
  igraph_empty(&h, 0, 0);         /* empty graph */
  igraph_vector_init(&quant, 1);  /* quantities vector */
  igraph_vector_init(&strat, 2);  /* strategies vector */

  /* test parameters */
  /*--------graph--vertex--optimality--quantities--strategies--mode------*/
  /* null pointer for graph */
  strategy_test_t null_graph = {NULL, 0, 0, NULL, NULL, IGRAPH_ALL};
  /* null pointer for quantities vector */
  strategy_test_t null_quant = {&g, 0, 0, NULL, NULL, IGRAPH_ALL};
  /* null pointer for strategies vector */
  strategy_test_t null_strat = {&g, 0, 0, &quant, NULL, IGRAPH_ALL};
  /* empty graph */
  strategy_test_t empty_graph = {&h, 0, 0, &quant, &strat, IGRAPH_ALL};
  /* length of quantities vector different from number of vertices */
  strategy_test_t qdiff_length = {&g, 0, 0, &quant, &strat, IGRAPH_ALL};
  /* length of strategies vector different from number of vertices */
  strategy_test_t sdiff_length = {&g, 0, 0, &quant, &strat, IGRAPH_ALL};
  strategy_test_t *all_checks[] = {/* 1 */ &null_graph,
				   /* 2 */ &null_quant,
				   /* 3 */ &null_strat,
				   /* 4 */ &empty_graph,
				   /* 5 */ &qdiff_length,
				   /* 6 */ &sdiff_length};
  n = 6;
  /* Run the error tests. We expect an error to be raised for each test. */
  igraph_set_error_handler(igraph_error_handler_ignore);
  i = 0;
  while (i < n) {
    test = all_checks[i];
    ret = igraph_deterministic_optimal_imitation(test->graph,
  						 test->vertex,
  						 test->optimality,
  						 test->quantities,
  						 test->strategies,
  						 test->mode);
    if (ret) {
      printf("Error test no. %d failed.\n", (int)(i + 1));
      return IGRAPH_FAILURE;
    }
    i++;
  }
  /* clean up */
  igraph_destroy(&g);
  igraph_destroy(&h);
  igraph_vector_destroy(&quant);
  igraph_vector_destroy(&strat);

  return IGRAPH_SUCCESS;
}

/* Updating the strategy of an isolated vertex. In this case, the strategies
 * vector should not change at all.
 */
int isolated_vertex_test() {
  igraph_t g;
  igraph_vector_t quant, strat, v;
  int i, ret;

  /* graph with one isolated vertex */
  igraph_vector_init(&v, 6);
  VECTOR(v)[0] = 0; VECTOR(v)[1] = 1;
  VECTOR(v)[2] = 1; VECTOR(v)[3] = 2;
  VECTOR(v)[4] = 2; VECTOR(v)[5] = 0;
  igraph_create(/*graph*/ &g, /*edge list*/ &v,
  		/*n vertices*/ 0, /*directed*/ 0);
  igraph_vector_destroy(&v);
  igraph_add_vertices(&g, 1, 0);  /* new vertex 3 is isolated */
  /* quantities vector: all vertices have the same fitness */
  igraph_vector_init(&quant, 4);
  VECTOR(quant)[0] = 0.25; VECTOR(quant)[1] = 0.25;
  VECTOR(quant)[2] = 0.25; VECTOR(quant)[3] = 0.25;
  /* strategies vector: 0 means aggressive strategy; 1 means passive */
  igraph_vector_init(&strat, 4);
  VECTOR(strat)[0] = 1; VECTOR(strat)[1] = 0;
  VECTOR(strat)[2] = 1; VECTOR(strat)[3] = 0;
  /* make a copy of the original strategies vector for comparison later on */
  igraph_vector_copy(&v, &strat);
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
  for (i = 0; i < igraph_vector_size(&strat); i++) {
    if (VECTOR(strat)[i] != VECTOR(v)[i]) {
      printf("Isolated vertex test failed.\n");
      return IGRAPH_FAILURE;
    }
  }
  /* clean up */
  igraph_destroy(&g);
  igraph_vector_destroy(&quant);
  igraph_vector_destroy(&strat);
  igraph_vector_destroy(&v);

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
int petersen_game_test() {
  igraph_t g;
  igraph_vector_t known_max_v, known_min_v, quant, strat, stratcopy, v;
  int i, nedge, nvert;

  /* the Petersen graph */
  nedge = 15;
  nvert = 10;
  igraph_vector_init(&v, 2 * nedge);
  VECTOR(v)[0]  = 0;  VECTOR(v)[1]  = 1;
  VECTOR(v)[2]  = 0;  VECTOR(v)[3]  = 4;
  VECTOR(v)[4]  = 0;  VECTOR(v)[5]  = 5;
  VECTOR(v)[6]  = 1;  VECTOR(v)[7]  = 2;
  VECTOR(v)[8]  = 1;  VECTOR(v)[9]  = 6;
  VECTOR(v)[10] = 2;  VECTOR(v)[11] = 3;
  VECTOR(v)[12] = 2;  VECTOR(v)[13] = 7;
  VECTOR(v)[14] = 3;  VECTOR(v)[15] = 4;
  VECTOR(v)[16] = 3;  VECTOR(v)[17] = 8;
  VECTOR(v)[18] = 4;  VECTOR(v)[19] = 9;
  VECTOR(v)[20] = 5;  VECTOR(v)[21] = 7;
  VECTOR(v)[22] = 5;  VECTOR(v)[23] = 8;
  VECTOR(v)[24] = 6;  VECTOR(v)[25] = 8;
  VECTOR(v)[26] = 6;  VECTOR(v)[27] = 9;
  VECTOR(v)[28] = 7;  VECTOR(v)[29] = 9;
  igraph_create(/*graph*/ &g, /*edge list*/ &v,
  		/*n vertices*/ 0, /*directed*/ 0);
  /* Strategies vector, one strategy for each vertex. Thus vec[i] is the */
  /* strategy of vertex i. The strategy space is: {0, 1, 2, 3}. */
  igraph_vector_init(&strat, nvert);
  VECTOR(strat)[0] = 1; VECTOR(strat)[1] = 1;
  VECTOR(strat)[2] = 2; VECTOR(strat)[3] = 2;
  VECTOR(strat)[4] = 0; VECTOR(strat)[5] = 0;
  VECTOR(strat)[6] = 0; VECTOR(strat)[7] = 1;
  VECTOR(strat)[8] = 2; VECTOR(strat)[9] = 3;
  /* Quantities vector, one quantity per vertex. Thus vec[i] is the */
  /* quantity for vertex i. */
  igraph_vector_init(&quant, nvert);
  VECTOR(quant)[0] = 0.3; VECTOR(quant)[1] = 1.1;
  VECTOR(quant)[2] = 0.5; VECTOR(quant)[3] = 1.0;
  VECTOR(quant)[4] = 0.9; VECTOR(quant)[5] = 0.8;
  VECTOR(quant)[6] = 0.4; VECTOR(quant)[7] = 0.1;
  VECTOR(quant)[8] = 0.7; VECTOR(quant)[9] = 0.7;
  /* Known strategies that would be adopted. Thus vec[i] means that in */
  /* game i where we revise the strategy of vertex i, the strategy */
  /* vec[i] would be adopted by i. */
  igraph_vector_init(&known_max_v, nvert);  /*maximum deterministic imitation*/
  VECTOR(known_max_v)[0] = 1; VECTOR(known_max_v)[1] = 1;
  VECTOR(known_max_v)[2] = 1; VECTOR(known_max_v)[3] = 2;
  VECTOR(known_max_v)[4] = 2; VECTOR(known_max_v)[5] = 0;
  VECTOR(known_max_v)[6] = 1; VECTOR(known_max_v)[7] = 0;
  VECTOR(known_max_v)[8] = 2; VECTOR(known_max_v)[9] = 0;
  igraph_vector_init(&known_min_v, nvert);  /*minimum deterministic imitation*/
  VECTOR(known_min_v)[0] = 1; VECTOR(known_min_v)[1] = 1;
  VECTOR(known_min_v)[2] = 1; VECTOR(known_min_v)[3] = 2;
  VECTOR(known_min_v)[4] = 1; VECTOR(known_min_v)[5] = 1;
  VECTOR(known_min_v)[6] = 0; VECTOR(known_min_v)[7] = 1;
  VECTOR(known_min_v)[8] = 0; VECTOR(known_min_v)[9] = 1;
  /* play game and compare resulting updated strategies */
  for (i = 0; i < nvert; i++) {
    /* maximum deterministic imitation */
    igraph_vector_copy(&stratcopy, &strat);
    igraph_deterministic_optimal_imitation(/*graph*/ &g,
					   /*vertex*/ (igraph_integer_t)i,
					   /*optimality*/ IGRAPH_MAXIMUM,
					   /*quantities*/ &quant,
					   /*strategies*/ &stratcopy,
					   /*neighbours*/ IGRAPH_ALL);
    if (VECTOR(stratcopy)[i] != VECTOR(known_max_v)[i]) {
      printf("Maximum deterministic imitation failed for vertex %d.\n", i);
      return IGRAPH_FAILURE;
    }
    igraph_vector_destroy(&stratcopy);
    /* minimum deterministic imitation */
    igraph_vector_copy(&stratcopy, &strat);
    igraph_deterministic_optimal_imitation(/*graph*/ &g,
					   /*vertex*/ (igraph_integer_t)i,
					   /*optimality*/ IGRAPH_MINIMUM,
					   /*quantities*/ &quant,
					   /*strategies*/ &stratcopy,
					   /*neighbours*/ IGRAPH_ALL);
    if (VECTOR(stratcopy)[i] != VECTOR(known_min_v)[i]) {
      printf("Minimum deterministic imitation failed for vertex %d.\n", i);
      return IGRAPH_FAILURE;
    }
    igraph_vector_destroy(&stratcopy);
  }
  /* clean up */
  igraph_destroy(&g);
  igraph_vector_destroy(&known_max_v);
  igraph_vector_destroy(&known_min_v);
  igraph_vector_destroy(&quant);
  igraph_vector_destroy(&strat);
  igraph_vector_destroy(&v);

  return IGRAPH_SUCCESS;
}

int main() {
  int ret;

  ret = error_tests();
  if (ret)
    return ret;
  ret = isolated_vertex_test();
  if (ret)
    return ret;
  ret = petersen_game_test();
  if (ret)
    return ret;

  return IGRAPH_SUCCESS;
}
