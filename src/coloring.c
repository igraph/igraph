
#include "igraph_coloring.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_interrupt_internal.h"
#include <assert.h>


int igraph_i_vertex_coloring_greedy_cn(const igraph_t *graph, igraph_vector_int_t *colors) {
    igraph_integer_t i, vertex, maxdeg;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_vector_int_t cn; /* number of already coloured neighbours */
    igraph_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_vector_int_resize(colors, vc));
    igraph_vector_int_fill(colors, 0);

    IGRAPH_CHECK(igraph_vector_int_init(&cn, vc));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &cn);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* find maximum degree and a corresponding vertex */
    {
        igraph_vector_t degree;

        IGRAPH_CHECK(igraph_vector_init(&degree, 0));
        IGRAPH_FINALLY(igraph_vector_destroy, &degree);
        IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, 0));

        vertex = (igraph_integer_t) igraph_vector_which_max(&degree);
        maxdeg = VECTOR(degree)[vertex];

        igraph_vector_destroy(&degree);
        IGRAPH_FINALLY_CLEAN(1);
    }

    while (VECTOR(cn)[vertex] >= 0) {
        igraph_vector_int_t *neighbors = igraph_adjlist_get(&adjlist, vertex);
        long neigh_count = igraph_vector_int_size(neighbors);

        /* colour current vertex */
        {
            igraph_integer_t col;
            igraph_vector_int_t neigh_colors;

            IGRAPH_CHECK(igraph_vector_int_init(&neigh_colors, neigh_count));
            for (i=0; i < neigh_count; ++i)
                VECTOR(neigh_colors)[i] = VECTOR(*colors)[ VECTOR(*neighbors)[i] ];
            igraph_vector_int_sort(&neigh_colors);

            i=0;
            col = 0;
            do {
                while (i < neigh_count && VECTOR(neigh_colors)[i] == col)
                    i++;
                col++;
            } while (i < neigh_count && VECTOR(neigh_colors)[i] == col);

            VECTOR(*colors)[vertex] = col;

            igraph_vector_int_destroy(&neigh_colors);
        }

        /* increment cn for all neighbours */
        for (i=0; i < neigh_count; ++i)
            VECTOR(cn)[ VECTOR(*neighbors)[i] ] += 1;

        /* assign a cn value to current vertex that will stay negative throughout the calculation */
        VECTOR(cn)[vertex] = -(maxdeg+1);

        /* choose next vertex with most already coloured neighbours */
        vertex = igraph_vector_int_which_max(&cn);

        IGRAPH_ALLOW_INTERRUPTION();
    }

    /* free data structures */
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&cn);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_vertex_coloring_greedy
 * \brief Computes a vertex colouring using a greedy algorithm.
 *
 * </para><para>
 * This function assigns a "colour"---represented as a positive integer---to
 * each vertex of the graph in such a way that neighbouring vertices never have
 * the same colour. The obtained colouring is not necessarily minimal.
 *
 * \param graph The input graph.
 * \param colors Pointer to an initialized integer vector. The colours will be stored here.
 * \param heuristic The vertex ordering heuristic to use during greedy colouring.
 *
 */
int igraph_vertex_coloring_greedy(const igraph_t *graph, igraph_vector_int_t *colors, igraph_coloring_greedy_t heuristic) {
    switch (heuristic) {
    case IGRAPH_COLORING_GREEDY_CN:
        return igraph_i_vertex_coloring_greedy_cn(graph, colors);
    }
}
