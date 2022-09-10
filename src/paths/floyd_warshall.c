
#include "igraph_paths.h"
#include "igraph_interface.h"

/**
 * \function igraph_distances_floyd_warshall
 * \brief Weighted all-pairs shortest path lengths with the Floyd-Warshall algorithm.
 *
 * The Floyd-Warshall algorithm computes weighted shortest path lengths between
 * all pairs of vertices at the same time. It is useful with very dense weighted graphs,
 * as its running time is primarily determined by the vertex count, and is not sensitive
 * to the graph density. In sparse graphs, other methods such as the Dijkstra or
 * Bellman-Ford algorithms will perform significantly better.
 *
 * \param graph The graph object.
 * \param res An intialized matrix, the distances will be stored here.
 * \param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \param weights The edge weights. If \c NULL, all weights are assumed to be 1.
 *   Negative weights are allowed, but the graph must not contain negative cycles.
 * \return Error code.
 *
 * \sa \ref igraph_distances(), \ref igraph_distances_dijkstra(),
 * \ref igraph_distances_bellman_ford(), \ref igraph_distances_johnson()
 *
 * Time complexity: O(|V|^3 + |E|) where |V| is the number of vertices
 * and |E| is the number of edges.
 */
igraph_error_t igraph_distances_floyd_warshall(
        const igraph_t *graph, igraph_matrix_t *res,
        igraph_neimode_t mode, const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t in = false, out = false;

    if (weights && igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERROR("Invalid weight vector length.", IGRAPH_EINVAL);
    }

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    switch (mode) {
    case IGRAPH_ALL:
        in = out = true;
        break;
    case IGRAPH_OUT:
        out = true;
        break;
    case IGRAPH_IN:
        in = true;
        break;
    default:
        IGRAPH_ERROR("Invalid mode.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_fill(res, IGRAPH_INFINITY);

    for (igraph_integer_t v=0; v < no_of_nodes; v++) {
        MATRIX(*res, v, v) = 0;
    }

    for (igraph_integer_t e=0; e < no_of_edges; e++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e);
        igraph_integer_t to = IGRAPH_TO(graph, e);
        igraph_real_t w = weights ? VECTOR(*weights)[e] : 1;

        if (weights && igraph_is_nan(w)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }

        if (w < 0) {
            if (mode == IGRAPH_ALL) {
                IGRAPH_ERRORF("Negative edge weight (%g) found in undirected graph "
                              "while calculating distances with Floyd-Warshall.",
                              IGRAPH_ENEGLOOP, w);
            } else if (to == from) {
                IGRAPH_ERRORF("Self-loop with negative weight (%g) found "
                              "while calculating distances with Floyd-Warshall.",
                              IGRAPH_ENEGLOOP, w);
            }
        }

        if (out && MATRIX(*res, from, to) > w) MATRIX(*res, from, to) = w;
        if (in  && MATRIX(*res, to, from) > w) MATRIX(*res, to, from) = w;
    }

    for (igraph_integer_t k=0; k < no_of_nodes; k++) {
        /* Iteration order matters for performance!
         * First j, then i, because matrices are stored as column-major. */
        for (igraph_integer_t j=0; j < no_of_nodes; j++) {
            for (igraph_integer_t i=0; i < no_of_nodes; i++) {
                igraph_real_t di = MATRIX(*res, i, k) + MATRIX(*res, k, j);
                igraph_real_t dd = MATRIX(*res, i, j);
                if (di < dd) {
                    MATRIX(*res, i, j) = di;
                }
                if (i == j && MATRIX(*res, i, i) < 0) {
                    IGRAPH_ERROR("Negative cycle found while calculating distances with Floyd-Warshall.",
                                 IGRAPH_ENEGLOOP);
                }
            }
        }
    }

    return IGRAPH_SUCCESS;
}
