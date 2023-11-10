#include "igraph_dispersion.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_dqueue.h"

#include "core/indheap.h"
#include "core/interruption.h"


igraph_error_t igraph_dispersion(const igraph_t *graph_u, igraph_integer_t u, igraph_integer_t v) {
    // Get the neighbors of u
    igraph_vector_int_t u_nbrs;
    igraph_vector_int_init(&u_nbrs, 0);
    igraph_neighbors(graph_u, &u_nbrs, u, IGRAPH_ALL);

    // Get the neighbors of v
    igraph_vector_int_t v_nbrs;
    igraph_vector_int_init(&v_nbrs, 0);
    igraph_neighbors(graph_u, &v_nbrs, v, IGRAPH_ALL);

    // Get the common neighbors of u and v
    igraph_vector_int_t u_v_nbrs;
    igraph_vector_int_init(&u_v_nbrs, 0);
    intersection(u_nbrs, v_nbrs, &u_v_nbrs);

    // Traverse all possible ties of connections that u and b share
    igraph_integer_t i, j, k, s, t, neighbors_size = igraph_vector_int_size(&u_v_nbrs);
    igraph_integer_t dispersion_val = 0;
    for (i = 0; i < neighbors_size - 1; i++) {
        // Get neighbor of s that are in G_u
        s = VECTOR(u_v_nbrs)[i];
        igraph_vector_int_t s_nbrs;
        igraph_vector_int_init(&s_nbrs, 0);
        igraph_neighbors(graph_u, &s_nbrs, s, IGRAPH_ALL);
        igraph_vector_int_t u_s_nbrs;
        igraph_vector_int_init(&u_s_nbrs, 0);
        intersection(u_nbrs, s_nbrs, &u_s_nbrs);
        for (j = i + 1; j < neighbors_size; j++) {
            t = VECTOR(u_v_nbrs)[j];
            // Determine if s and t are directly connected
            if (igraph_vector_int_contains(&u_s_nbrs, t)) {
                continue;
            }
            // Determine if the path from s to t is 2
            igraph_vector_int_t t_nbrs;
            igraph_vector_int_init(&t_nbrs, 0);
            igraph_neighbors(graph_u, &t_nbrs, t, IGRAPH_ALL);
            // Determine if s and t share a connection
            for (k = 0; k < igraph_vector_int_size(&u_s_nbrs); k++) {
                // Skip u and v
                if (VECTOR(u_s_nbrs)[k] == u || VECTOR(u_s_nbrs)[k] == v) continue;
                if (igraph_vector_int_contains(&t_nbrs, VECTOR(u_s_nbrs)[k])) {
                    break;
                }
            }
            if (k == igraph_vector_int_size(&u_s_nbrs)) dispersion_val++;
        }
    }
    return dispersion_val;
}

igraph_error_t intersection(igraph_vector_int_t u_nbrs, igraph_vector_int_t v_nbrs, igraph_vector_int_t *u_v_nbrs) {
    igraph_integer_t u_nbrs_num, v_nbrs_num, i, j, pos = 0;
    u_nbrs_num = igraph_vector_int_size(&u_nbrs);
    v_nbrs_num = igraph_vector_int_size(&v_nbrs);
    IGRAPH_CHECK(igraph_vector_int_resize(u_v_nbrs, u_nbrs_num));
    for (i = 0; i < u_nbrs_num; i++) {
        for (j = 0; j < v_nbrs_num; j++) {
            if (VECTOR(u_nbrs)[i] == VECTOR(v_nbrs)[j]) {
                VECTOR(*u_v_nbrs)[pos++] = VECTOR(u_nbrs)[i];
            }
        }
    }
    IGRAPH_CHECK(igraph_vector_int_resize(u_v_nbrs, pos));
    return 0;
}


