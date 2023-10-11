//
// Created by Lara Holm on 11.10.2023.
//
//
// Created by Lara Holm on 18.9.2023.
//
#include <igraph.h>
#include <stdio.h>


/**
 * Populates a Joint Degree Matrix. Works on directed and undirected graphs.
 * @param g_ptr A pointer to an initialized graph object.
 * @param jdm_ptr A pointer to a zero initialized integer matrix. The function writes the appropriate values here.
 * @param degrees_ptr A pointer to an initialized integer vector containing the degrees of the graph object.
 * @return 0
 */
int construct_jdm(igraph_t* g_ptr, igraph_matrix_int_t* jdm_ptr, igraph_vector_int_t* degrees_ptr) {
    // Get the edges
    igraph_eit_t eit;
    igraph_es_t es;
    igraph_integer_t eid;
    igraph_integer_t v1id;
    igraph_integer_t v2id;
    igraph_integer_t v1deg;
    igraph_integer_t v2deg;

    igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
    igraph_eit_create(g_ptr, es, &eit);

    if (igraph_is_directed(g_ptr)) {
        printf("Directed\n");

        while (!IGRAPH_EIT_END(eit)) {
            v1id = IGRAPH_FROM(g_ptr, IGRAPH_EIT_GET(eit));
            v2id = IGRAPH_TO(g_ptr, IGRAPH_EIT_GET(eit));
            eid = IGRAPH_EIT_GET(eit);
            igraph_degree_1(g_ptr, &v1deg, v1id, IGRAPH_OUT, false);
            igraph_degree_1(g_ptr, &v2deg, v2id, IGRAPH_IN, false);

            // If the graph is directed, J out_degree, in_degree
            // MATRIX(m , i, j) = macro in igraph
            MATRIX(*jdm_ptr, v1deg, v2deg) = MATRIX(*jdm_ptr, v1deg, v2deg) + 1;
            // igraph_matrix_int_set(jdm_ptr, v1deg, v2deg, igraph_matrix_int_get(jdm_ptr, v1deg, v2deg) + 1);
            printf("EdgeID: %" IGRAPH_PRId, IGRAPH_EIT_GET(eit));
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
            v1id = IGRAPH_FROM(g_ptr, IGRAPH_EIT_GET(eit));
            v2id = IGRAPH_TO(g_ptr, IGRAPH_EIT_GET(eit));
            eid = IGRAPH_EIT_GET(eit);
            v1deg = igraph_vector_int_get(degrees_ptr, v1id);
            v2deg = igraph_vector_int_get(degrees_ptr, v2id);

            // If the graph is undirected, it is symmetrical and the sum of the entire matrix is 2x the number of edges.
            //igraph_matrix_int_set(jdm_ptr, v1deg, v2deg, igraph_matrix_int_get(jdm_ptr, v1deg, v2deg) + 1);
            //igraph_matrix_int_set(jdm_ptr, v2deg, v1deg, igraph_matrix_int_get(jdm_ptr, v2deg, v1deg) + 1);
            MATRIX(*jdm_ptr, v1deg, v2deg) = MATRIX(*jdm_ptr, v1deg, v2deg) + 1;
            MATRIX(*jdm_ptr, v2deg, v1deg) = MATRIX(*jdm_ptr, v2deg, v1deg) + 1;

            printf("EdgeID: %" IGRAPH_PRId, IGRAPH_EIT_GET(eit));
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

    return 0;
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
    igraph_matrix_int_init(&jdm, g_max_degree + 1, g_max_degree + 1);
    igraph_matrix_int_init(&jdm_dir, g_dir_max_out_degree+ 1, g_dir_max_in_degree + 1);


    construct_jdm(&g, &jdm, &g_degrees);
    construct_jdm(&g_dir, &jdm_dir, &g_dir_degrees);

    igraph_matrix_int_print(&jdm);
    igraph_matrix_int_print(&jdm_dir);

    // Clean up
    igraph_destroy(&g);
    igraph_destroy(&g_dir);
    igraph_matrix_int_destroy(&jdm);
    igraph_matrix_int_destroy(&jdm_dir);

    igraph_vector_int_destroy(&g_degrees);
    igraph_vector_int_destroy(&g_dir_degrees);

    return 0;
}
