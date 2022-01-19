#include <igraph.h>
#include "test_utilities.inc"

/* --------------------------------------------------------- */
/* new typedef proposed to add to include/igraph_constants.h */
/* --------------------------------------------------------- */

typedef enum { IGRAPH_WHEEL_OUT = 0, IGRAPH_WHEEL_IN,
               IGRAPH_WHEEL_UNDIRECTED,
               IGRAPH_WHEEL_MUTUAL
             } igraph_wheel_mode_t;

/* -------------------------------------------------------------- */
/* new igraph_wheel proposed to add to src/constructors/regualr.c */
/* -------------------------------------------------------------- */
/**
 * \ingroup generators
 * \function igraph_wheel
 * \brief Creates a \em wheel graph, a union of a star graph and a
 *        cycle graph.
 * \param graph Pointer to an uninitialized graph object, this will
 *        be the result.
 * \param n Integer constant, the number of vertices in the graph.
 * \param mode Constant, gives the type of the star graph to
 *        create. Possible values:
 *        \clist
 *        \cli IGRAPH_WHEEL_OUT
 *          directed wheel graph, edges point
 *          \em from the center to the other vertices.
 *        \cli IGRAPH_WHEEL_IN
 *          directed wheel graph, edges point
 *          \em to the center from the other vertices.
 *        \cli IGRAPH_STAR_MUTUAL
 *          directed wheel graph with mutual edges.
 *        \cli IGRAPH_STAR_UNDIRECTED
 *          an undirected wheel graph is
 *          created.
 *        \endclist
 * \param center Id of the vertex which will be the center of the
 *          graph.
 * \return Error code:
 *         \clist
 *         \cli IGRAPH_EINVVID
 *           invalid number of vertices.
 *         \cli IGRAPH_EINVAL
 *           invalid center vertex.
 *         \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *         \endclist
 *
 * Time complexity: O(|V|), the
 * number of vertices in the graph.
 *
 * \sa \ref igraph_lattice(), \ref igraph_ring(), \ref igraph_star, \ref igraph_tree()
 * for creating other regular structures.
 *
 */

int igraph_wheel(igraph_t *graph, igraph_integer_t n, igraph_wheel_mode_t mode,
                igraph_integer_t center) {

    igraph_integer_t* rim_vertex;
    igraph_integer_t  v;
    igraph_star_mode_t star_mode;
    igraph_vector_t rim_edges;
    long int i;

    /* create a star and also make use of its existing pre-check */
    /* works for all input parameters                            */

    switch(mode)
    {
        case IGRAPH_WHEEL_OUT:
            star_mode = IGRAPH_STAR_OUT;
            break;
        case IGRAPH_WHEEL_IN:
            star_mode = IGRAPH_STAR_IN;
            break;
        case IGRAPH_WHEEL_MUTUAL:
            star_mode = IGRAPH_STAR_MUTUAL;
            break;
        case IGRAPH_WHEEL_UNDIRECTED:
            star_mode = IGRAPH_STAR_UNDIRECTED;
            break;
        default:
            IGRAPH_ERROR("invalid mode", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_star(graph, n, star_mode, center));

    /* If n <= 2, wheel graph is identical with star graph, */
    /* no further processing is needed                      */

    if (n <= 2) {
        return 0;
    }

    /* A wheel rim will be created with an array rim_vertex */
    /* whose elements contains the real vertices            */

    rim_vertex = (igraph_integer_t*)calloc(n-2, sizeof(igraph_integer_t));
    for (i = 0; i < center; i++) {
        rim_vertex[i] = i;
    }
    for (i = center; i < n-1; i++) {
        rim_vertex[i] = i + 1;
    }

    /* Add edges to the rim */

    igraph_vector_init(&rim_edges, 2 * (n-1));
    for (i = 0; i < n-2; i++) {
        VECTOR(rim_edges)[2 * i] = rim_vertex[i];
        VECTOR(rim_edges)[2 * i + 1] = rim_vertex[i + 1];
    }
    VECTOR(rim_edges)[2 * n - 4] = rim_vertex[n-2];
    VECTOR(rim_edges)[2 * n - 3] = rim_vertex[0];

    /* Combine the rim into the star to make it a wheel graph */

    IGRAPH_CHECK(igraph_add_edges(graph, &rim_edges, 0));

    /* Add a reverse direction rim if mode is MUTUAL */

    if (mode == IGRAPH_WHEEL_MUTUAL) {
        for (i=0; i < n-1; i++) {
            v = VECTOR(rim_edges)[2 * i];
            VECTOR(rim_edges)[2 * i] = VECTOR(rim_edges)[2 * i + 1];
            VECTOR(rim_edges)[2 * i + 1] = v;
        }
        IGRAPH_CHECK(igraph_add_edges(graph, &rim_edges, 0));
    }

    igraph_vector_destroy(&rim_edges);
    free(rim_vertex);

    return 0;
    
}

/* -------------------------------------------------------------- */
/* unit test of igraph_wheel */
/* -------------------------------------------------------------- */

void call_and_print(igraph_t *graph, igraph_integer_t n, igraph_wheel_mode_t mode,
                igraph_integer_t center) {

    IGRAPH_ASSERT(igraph_wheel(graph, n, mode, center) == IGRAPH_SUCCESS);
    igraph_write_graph_edgelist(graph, stdout);
    igraph_destroy(graph);
    printf("\n");
}

int main() {
    igraph_t graph;

    printf("1 vertex:\n");
    call_and_print(&graph, 1, IGRAPH_STAR_UNDIRECTED, 0);
    printf("2 vertices:\n");
    call_and_print(&graph, 2, IGRAPH_STAR_UNDIRECTED, 0);
    printf("OUT:\n");
    call_and_print(&graph, 4, IGRAPH_STAR_OUT, 0);
    printf("IN:\n");
    call_and_print(&graph, 4, IGRAPH_STAR_IN, 0);
    printf("MUTUAL:\n");
    call_and_print(&graph, 4, IGRAPH_STAR_MUTUAL, 0);
    printf("UNDIRECTED:\n");
    call_and_print(&graph, 4, IGRAPH_STAR_UNDIRECTED, 0);
    printf("center is n/2:\n");
    call_and_print(&graph, 4, IGRAPH_STAR_OUT, 2);
    printf("center is n - 1:\n");
    call_and_print(&graph, 4, IGRAPH_STAR_OUT, 3);

    return 0;
}
