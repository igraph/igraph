/*
   IGraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>
#include "test_utilities.h"


/**
 * \function igraph_construct_jdm
 * \brief Constructs a joint degree matrix.
 *
 * The joint degree matrix of a graph contains the number of edges between vertices of degree i and degree j for every
 * (i, j). The function populates a joint degree matrix and works on directed and undirected graphs.
 *
 * \param graph A pointer to an initialized graph object.
 * \param jdm_ptr A pointer to a zero initialized integer matrix. The function writes the appropriate values here.
 * \param degrees_ptr A pointer to an initialized integer vector containing the degrees of the graph object.
 * \return Error code.
 *
 * TODO: Time complexity: O(V^2 d),
 * where V is the number of vertices in the vertex iterator given, and d is the
 * (maximum) degree of the vertices in the graph.
 *
 * \sa \ref igraph_similarity_jaccard(), a measure very similar to the Dice
 *   coefficient.
 *
 * \example examples/simple/igraph_similarity.c
 */
igraph_error_t igraph_construct_jdm(igraph_t* graph, igraph_matrix_int_t* jdm_ptr, igraph_vector_int_t* degrees_ptr) {
    // Get the edges
    igraph_eit_t eit;
    igraph_es_t es;
    igraph_integer_t eid;
    igraph_integer_t v1id;
    igraph_integer_t v2id;
    igraph_integer_t v1deg;
    igraph_integer_t v2deg;
    igraph_vector_int_t degrees;
    igraph_vector_int_t out_degrees;
    igraph_vector_int_t in_degrees;

    IGRAPH_CHECK(igraph_vector_int_init(&degrees, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&out_degrees, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&in_degrees, 0));

    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, false));
    IGRAPH_CHECK(igraph_degree(graph, &out_degrees, igraph_vss_all(), IGRAPH_OUT, false));
    IGRAPH_CHECK(igraph_degree(graph, &in_degrees, igraph_vss_all(), IGRAPH_IN, false));

    IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));

    if (igraph_is_directed(graph)) {
        printf("Directed\n");

        while (!IGRAPH_EIT_END(eit)) {
            eid = IGRAPH_EIT_GET(eit);
            v1id = IGRAPH_FROM(graph, eid);
            v2id = IGRAPH_TO(graph, eid);
            //igraph_degree_1(graph, &v1deg, v1id, IGRAPH_OUT, false);
            //igraph_degree_1(graph, &v2deg, v2id, IGRAPH_IN, false);
            //
            v1deg = igraph_vector_int_get(&out_degrees,v1id);
            v2deg = igraph_vector_int_get(&in_degrees, v2id);

            // If the graph is directed, J out_degree, in_degree
            // MATRIX(m , i, j) = macro in igraph
            MATRIX(*jdm_ptr, v1deg-1, v2deg-1)++; //= MATRIX(*jdm_ptr, v1deg, v2deg) + 1;
            // igraph_matrix_int_set(jdm_ptr, v1deg, v2deg, igraph_matrix_int_get(jdm_ptr, v1deg, v2deg) + 1);
            printf("EdgeID: %" IGRAPH_PRId, eid);
            printf("\n");
            printf(" %" IGRAPH_PRId, v1id);
            printf(" Degree: %" IGRAPH_PRId, v1deg);
            printf("\n");
            printf(" %" IGRAPH_PRId, v2id);
            printf(" Degree: %" IGRAPH_PRId, v2deg);
            printf("\n");
            IGRAPH_EIT_NEXT(eit);
        }
        printf("\n");
    } else {
        printf("Undirected\n");

        while (!IGRAPH_EIT_END(eit)) {
            eid = IGRAPH_EIT_GET(eit);
            v1id = IGRAPH_FROM(graph, eid);
            v2id = IGRAPH_TO(graph, eid);

//            v1deg = igraph_vector_int_get(degrees_ptr, v1id);
//            v2deg = igraph_vector_int_get(degrees_ptr, v2id);

            v1deg = igraph_vector_int_get(&degrees, v1id);
            v2deg = igraph_vector_int_get(&degrees, v2id);

            // If the graph is undirected, it is symmetrical and the sum of the entire matrix is 2x the number of edges.
            //igraph_matrix_int_set(jdm_ptr, v1deg, v2deg, igraph_matrix_int_get(jdm_ptr, v1deg, v2deg) + 1);
            //igraph_matrix_int_set(jdm_ptr, v2deg, v1deg, igraph_matrix_int_get(jdm_ptr, v2deg, v1deg) + 1);
            MATRIX(*jdm_ptr, v1deg-1, v2deg-1)++; //= MATRIX(*jdm_ptr, v1deg, v2deg) + 1;
            MATRIX(*jdm_ptr, v2deg-1, v1deg-1)++; //= MATRIX(*jdm_ptr, v2deg, v1deg) + 1;

            printf("EdgeID: %" IGRAPH_PRId, eid);
            printf("\n");
            printf(" %" IGRAPH_PRId, v1id);
            printf(" Degree: %" IGRAPH_PRId, v1deg);
            printf("\n");
            printf(" %" IGRAPH_PRId, v2id);
            printf(" Degree: %" IGRAPH_PRId, v2deg);
            printf("\n");
            IGRAPH_EIT_NEXT(eit);
        }
        printf("\n");
    }
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&out_degrees);
    igraph_vector_int_destroy(&in_degrees);


    return IGRAPH_SUCCESS;
}

int main (void) {
    // Structures that need to be destroyed
    igraph_t g;
    igraph_t g_dir;
    igraph_matrix_int_t jdm;
    igraph_matrix_int_t jdm_dir;
    igraph_vector_int_t g_degrees;
    igraph_vector_int_t g_dir_degrees;

    igraph_integer_t g_max_degree;
    igraph_integer_t g_dir_max_out_degree;
    igraph_integer_t g_dir_max_in_degree;
    igraph_bool_t is_directed;

    igraph_small(&g, 5, false,
                 0, 1, 0, 2, 0, 4,
                 1, 4,
                 2, 3, 2, 4,
                 3, 4,
                 -1);

    igraph_small(&g_dir, 5, true,
                 0, 1, 0, 2, 0, 4,
                 1, 0,
                 3, 2, 3, 4,
                 4, 0, 4, 1, 4, 2,
                 -1);

//    igraph_rng_seed(igraph_rng_default(), 42);
//
//    igraph_erdos_renyi_game_gnm(
//            &g, 5, 7,
//            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS
//    );
//
//    igraph_rng_seed(igraph_rng_default(), 42);
//
//    igraph_erdos_renyi_game_gnm(
//            &g_dir, 5, 7,
//            IGRAPH_DIRECTED, IGRAPH_NO_LOOPS
//    );

    printf("===UNDIRECTED===:\n");
    printf("vcount: %d, ecount: %d\n\n", (int) igraph_vcount(&g), (int) igraph_ecount(&g));

    printf("===DIRECTED===:\n");
    printf("vcount: %d, ecount: %d\n\n", (int) igraph_vcount(&g_dir), (int) igraph_ecount(&g_dir));

    is_directed = igraph_is_directed(&g);
    printf("Is directed: %s\n", is_directed ? "true" : "false");

    is_directed = igraph_is_directed(&g_dir);
    printf("Is directed: %s\n", is_directed ? "true" : "false");

    // Compute and print the degrees:
    // UNDIRECTED
    igraph_vector_int_init(&g_degrees, 0);
    igraph_degree(&g, &g_degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
    // Crash if vector is empty - igraph_vector_* functions
    igraph_maxdegree(&g, &g_max_degree, igraph_vss_all(), IGRAPH_ALL, false);

    // DIRECTED
    igraph_vector_int_init(&g_dir_degrees, 0);

    igraph_degree(&g_dir, &g_dir_degrees, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS_TWICE);

    // Crash if vector is empty - igraph_vector_* functions
    igraph_maxdegree(&g_dir, &g_dir_max_in_degree, igraph_vss_all(), IGRAPH_IN, false);
    igraph_maxdegree(&g_dir, &g_dir_max_out_degree, igraph_vss_all(), IGRAPH_OUT, false);
    //printf("InDegrees: %d\n", (int)g_dir_max_in_degree);
    //printf("OutDegrees: %d\n", (int)g_dir_max_out_degree);


//    igraph_vector_int_print(&g_degrees);
//    printf("Max_degree: %d\n", (int) g_max_degree);

    // Zero initialize a matrix
    // Note: +1 if 0 row col should be included. If removed, remember to subtract 1 from the indices in the while loops.
    // initialize outside of the function, igraph_matrix_resize
    igraph_matrix_int_init(&jdm, g_max_degree, g_max_degree);
    igraph_matrix_int_init(&jdm_dir, g_dir_max_out_degree, g_dir_max_in_degree);


    igraph_construct_jdm(&g, &jdm, &g_degrees);
    igraph_construct_jdm(&g_dir, &jdm_dir, &g_dir_degrees);

    igraph_matrix_int_print(&jdm);
    igraph_matrix_int_print(&jdm_dir);

    // Clean up
    igraph_destroy(&g);
    igraph_destroy(&g_dir);
    igraph_matrix_int_destroy(&jdm);
    igraph_matrix_int_destroy(&jdm_dir);

    igraph_vector_int_destroy(&g_degrees);
    igraph_vector_int_destroy(&g_dir_degrees);

    VERIFY_FINALLY_STACK();

    return 0;
}
