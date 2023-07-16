/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_paths.h"

#include "igraph_conversion.h"
#include "igraph_interface.h"

#include "math/safe_intop.h"
/**
 * \function igraph_distances_johnson
 * \brief Weighted shortest path lengths between vertices, using Johnson's algorithm.
 *
 * This algorithm supports directed graphs with negative edge weights, and performs
 * better than the Bellman-Ford method when distances are calculated from many different
 * sources, the typical use case being all-pairs distance calculations. It works by using
 * a single-source Bellman-Ford run to transform all edge weights to non-negative ones,
 * then invoking Dijkstra's algorithm with the new weights. See the Wikipedia page
 * for more details: http://en.wikipedia.org/wiki/Johnson's_algorithm.
 *
 * </para><para>
 * If no edge weights are supplied, then the unweighted version, \ref igraph_distances()
 * is called. If none of the supplied edge weights are negative, then Dijkstra's algorithm
 * is used by calling \ref igraph_distances_dijkstra().
 *
 * </para><para>
 * Note that Johnson's algorithm applies only to directed graphs. This function rejects
 * undirected graphs with \em any negative edge weights, even when the \p from and \p to
 * vertices are all in connected components that are free of negative weights.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Donald B. Johnson: Efficient Algorithms for Shortest Paths in Sparse Networks.
 * J. ACM 24, 1 (1977), 1–13.
 * https://doi.org/10.1145/321992.321993
 *
 * \param graph The input graph. If negative weights are present, it
 *   should be directed.
 * \param res Pointer to an initialized matrix, the result will be
 *   stored here, one line for each source vertex, one column for each
 *   target vertex.
 * \param from The source vertices.
 * \param to The target vertices. It is not allowed to include a
 *   vertex twice or more.
 * \param weights Optional edge weights. If it is a null-pointer, then
 *   the unweighted breadth-first search based \ref igraph_distances() will
 *   be called. Edges with positive infinite weights are ignored.
 * \return Error code.
 *
 * Time complexity: O(s|V|log|V|+|V||E|), |V| and |E| are the number
 * of vertices and edges, s is the number of source vertices.
 *
 * \sa \ref igraph_distances() for a faster unweighted version,
 * \ref igraph_distances_dijkstra() if you do not have negative
 * edge weights, \ref igraph_distances_bellman_ford() if you only
 * need to calculate shortest paths from a couple of sources.
 */
igraph_error_t igraph_distances_johnson(const igraph_t *graph,
                                  igraph_matrix_t *res,
                                  const igraph_vs_t from,
                                  const igraph_vs_t to,
                                  const igraph_vector_t *weights) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_t newgraph;
    igraph_vector_int_t edges;
    igraph_vector_t newweights;
    igraph_matrix_t bfres;
    igraph_integer_t i, ptr;
    igraph_integer_t nr, nc;
    igraph_vit_t fromvit;
    igraph_integer_t no_edges_reserved;

    /* If no weights, then we can just run the unweighted version */
    if (!weights) {
        return igraph_distances(graph, res, from, to, IGRAPH_OUT);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    /* If no edges, then we can just run the unweighted version */
    if (no_of_edges == 0) {
        return igraph_distances(graph, res, from, to, IGRAPH_OUT);
    }

    /* If no negative weights, then we can run Dijkstra's algorithm */
    {
        igraph_real_t min_weight = igraph_vector_min(weights);
        if (isnan(min_weight)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
        if (min_weight >= 0) {
            return igraph_distances_dijkstra(graph, res, from, to, weights, IGRAPH_OUT);
        }
    }

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Johnson's shortest path: undirected graph and negative weight.",
                     IGRAPH_EINVAL);
    }

    /* ------------------------------------------------------------ */
    /* -------------------- Otherwise proceed --------------------- */

    IGRAPH_MATRIX_INIT_FINALLY(&bfres, 0, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&newweights, 0);

    IGRAPH_CHECK(igraph_empty(&newgraph, no_of_nodes + 1, igraph_is_directed(graph)));
    IGRAPH_FINALLY(igraph_destroy, &newgraph);

    IGRAPH_SAFE_MULT(no_of_nodes, 2, &no_edges_reserved);
    IGRAPH_SAFE_ADD(no_edges_reserved, no_of_edges * 2, &no_edges_reserved);

    /* Add a new node to the graph, plus edges from it to all the others. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_edges_reserved);
    igraph_get_edgelist(graph, &edges, /*bycol=*/ 0); /* reserved */
    igraph_vector_int_resize(&edges, no_edges_reserved); /* reserved */
    for (i = 0, ptr = no_of_edges * 2; i < no_of_nodes; i++) {
        VECTOR(edges)[ptr++] = no_of_nodes;
        VECTOR(edges)[ptr++] = i;
    }
    IGRAPH_CHECK(igraph_add_edges(&newgraph, &edges, 0));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_vector_reserve(&newweights, no_of_edges + no_of_nodes));
    igraph_vector_update(&newweights, weights); /* reserved */
    igraph_vector_resize(&newweights, no_of_edges + no_of_nodes); /* reserved */
    for (i = no_of_edges; i < no_of_edges + no_of_nodes; i++) {
        VECTOR(newweights)[i] = 0;
    }

    /* Run Bellman-Ford algorithm on the new graph, starting from the
       new vertex.  */

    IGRAPH_CHECK(igraph_distances_bellman_ford(&newgraph, &bfres,
                 igraph_vss_1(no_of_nodes),
                 igraph_vss_all(), &newweights, IGRAPH_OUT));

    igraph_destroy(&newgraph);
    IGRAPH_FINALLY_CLEAN(1);

    /* Now the edges of the original graph are reweighted, using the
       values from the BF algorithm. Instead of w(u,v) we will have
       w(u,v) + h(u) - h(v) */

    igraph_vector_resize(&newweights, no_of_edges); /* reserved */
    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t ffrom = IGRAPH_FROM(graph, i);
        igraph_integer_t tto = IGRAPH_TO(graph, i);
        VECTOR(newweights)[i] += MATRIX(bfres, 0, ffrom) - MATRIX(bfres, 0, tto);

        /* If a weight becomes slightly negative due to roundoff errors,
           snap it to exact zero. */
        if (VECTOR(newweights)[i] < 0) VECTOR(newweights)[i] = 0;
    }

    /* Run Dijkstra's algorithm on the new weights */
    IGRAPH_CHECK(igraph_distances_dijkstra(graph, res, from,
                 to, &newweights,
                 IGRAPH_OUT));

    igraph_vector_destroy(&newweights);
    IGRAPH_FINALLY_CLEAN(1);

    /* Reweight the shortest paths */
    nr = igraph_matrix_nrow(res);
    nc = igraph_matrix_ncol(res);

    IGRAPH_CHECK(igraph_vit_create(graph, from, &fromvit));
    IGRAPH_FINALLY(igraph_vit_destroy, &fromvit);

    for (i = 0; i < nr; i++, IGRAPH_VIT_NEXT(fromvit)) {
        igraph_integer_t v1 = IGRAPH_VIT_GET(fromvit);
        if (igraph_vs_is_all(&to)) {
            igraph_integer_t v2;
            for (v2 = 0; v2 < nc; v2++) {
                igraph_real_t sub = MATRIX(bfres, 0, v1) - MATRIX(bfres, 0, v2);
                MATRIX(*res, i, v2) -= sub;
            }
        } else {
            igraph_integer_t j;
            igraph_vit_t tovit;
            IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
            IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
            for (j = 0, IGRAPH_VIT_RESET(tovit); j < nc; j++, IGRAPH_VIT_NEXT(tovit)) {
                igraph_integer_t v2 = IGRAPH_VIT_GET(tovit);
                igraph_real_t sub = MATRIX(bfres, 0, v1) - MATRIX(bfres, 0, v2);
                MATRIX(*res, i, j) -= sub;
            }
            igraph_vit_destroy(&tovit);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    igraph_vit_destroy(&fromvit);
    igraph_matrix_destroy(&bfres);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_shortest_paths_johnson
 * \brief Weighted shortest path lengths between vertices, using Johnson's algorithm (deprecated).
 *
 * \deprecated-by igraph_distances_johnson 0.10.0
 */
igraph_error_t igraph_shortest_paths_johnson(const igraph_t *graph,
                                  igraph_matrix_t *res,
                                  const igraph_vs_t from,
                                  const igraph_vs_t to,
                                  const igraph_vector_t *weights) {
    return igraph_distances_johnson(graph, res, from, to, weights);
}
