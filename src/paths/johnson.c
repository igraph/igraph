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

/**
 * \function igraph_shortest_paths_johnson
 * \brief Weighted shortest path lengths between vertices, using Johnson's algorithm.
 *
 * See Wikipedia at http://en.wikipedia.org/wiki/Johnson's_algorithm
 * for Johnson's algorithm. This algorithm works even if the graph
 * contains negative edge weights, and it is worth using it if we
 * calculate the shortest paths from many sources.
 *
 * </para><para>
 * If no edge weights are supplied, then the unweighted
 * version, \ref igraph_shortest_paths() is called.
 *
 * </para><para>
 * If all the supplied edge weights are non-negative,
 * then Dijkstra's algorithm is used by calling
 * \ref igraph_shortest_paths_dijkstra().
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
 *   the unweighted breadth-first search based \ref
 *   igraph_shortest_paths() will be called.
 * \return Error code.
 *
 * Time complexity: O(s|V|log|V|+|V||E|), |V| and |E| are the number
 * of vertices and edges, s is the number of source vertices.
 *
 * \sa \ref igraph_shortest_paths() for a faster unweighted version
 * or \ref igraph_shortest_paths_dijkstra() if you do not have negative
 * edge weights, \ref igraph_shortest_paths_bellman_ford() if you only
 * need to calculate shortest paths from a couple of sources.
 */
int igraph_shortest_paths_johnson(const igraph_t *graph,
                                  igraph_matrix_t *res,
                                  const igraph_vs_t from,
                                  const igraph_vs_t to,
                                  const igraph_vector_t *weights) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_t newgraph;
    igraph_vector_t edges, newweights;
    igraph_matrix_t bfres;
    long int i, ptr;
    long int nr, nc;
    igraph_vit_t fromvit;

    /* If no weights, then we can just run the unweighted version */
    if (!weights) {
        return igraph_shortest_paths(graph, res, from, to, IGRAPH_OUT);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
    }

    /* If no edges, then we can just run the unweighted version */
    if (no_of_edges == 0) {
        return igraph_shortest_paths(graph, res, from, to, IGRAPH_OUT);
    }

    /* If no negative weights, then we can run Dijkstra's algorithm */
    {
        igraph_real_t min_weight = igraph_vector_min(weights);
        if (igraph_is_nan(min_weight)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values", IGRAPH_EINVAL);
        }
        if (min_weight >= 0) {
            return igraph_shortest_paths_dijkstra(graph, res, from, to,
                                                  weights, IGRAPH_OUT);
        }
    }

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("Johnson's shortest path: undirected graph and negative weight",
                     IGRAPH_EINVAL);
    }

    /* ------------------------------------------------------------ */
    /* -------------------- Otherwise proceed --------------------- */

    IGRAPH_MATRIX_INIT_FINALLY(&bfres, 0, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&newweights, 0);

    IGRAPH_CHECK(igraph_empty(&newgraph, (igraph_integer_t) no_of_nodes + 1,
                              igraph_is_directed(graph)));
    IGRAPH_FINALLY(igraph_destroy, &newgraph);

    /* Add a new node to the graph, plus edges from it to all the others. */
    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2 + no_of_nodes * 2);
    igraph_get_edgelist(graph, &edges, /*bycol=*/ 0);
    igraph_vector_resize(&edges, no_of_edges * 2 + no_of_nodes * 2);
    for (i = 0, ptr = no_of_edges * 2; i < no_of_nodes; i++) {
        VECTOR(edges)[ptr++] = no_of_nodes;
        VECTOR(edges)[ptr++] = i;
    }
    IGRAPH_CHECK(igraph_add_edges(&newgraph, &edges, 0));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_vector_reserve(&newweights, no_of_edges + no_of_nodes));
    igraph_vector_update(&newweights, weights);
    igraph_vector_resize(&newweights, no_of_edges + no_of_nodes);
    for (i = no_of_edges; i < no_of_edges + no_of_nodes; i++) {
        VECTOR(newweights)[i] = 0;
    }

    /* Run Bellmann-Ford algorithm on the new graph, starting from the
       new vertex.  */

    IGRAPH_CHECK(igraph_shortest_paths_bellman_ford(&newgraph, &bfres,
                 igraph_vss_1((igraph_integer_t) no_of_nodes),
                 igraph_vss_all(), &newweights, IGRAPH_OUT));

    igraph_destroy(&newgraph);
    IGRAPH_FINALLY_CLEAN(1);

    /* Now the edges of the original graph are reweighted, using the
       values from the BF algorithm. Instead of w(u,v) we will have
       w(u,v) + h(u) - h(v) */

    igraph_vector_resize(&newweights, no_of_edges);
    for (i = 0; i < no_of_edges; i++) {
        long int ffrom = IGRAPH_FROM(graph, i);
        long int tto = IGRAPH_TO(graph, i);
        VECTOR(newweights)[i] += MATRIX(bfres, 0, ffrom) - MATRIX(bfres, 0, tto);
    }

    /* Run Dijkstra's algorithm on the new weights */
    IGRAPH_CHECK(igraph_shortest_paths_dijkstra(graph, res, from,
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
        long int v1 = IGRAPH_VIT_GET(fromvit);
        if (igraph_vs_is_all(&to)) {
            long int v2;
            for (v2 = 0; v2 < nc; v2++) {
                igraph_real_t sub = MATRIX(bfres, 0, v1) - MATRIX(bfres, 0, v2);
                MATRIX(*res, i, v2) -= sub;
            }
        } else {
            long int j;
            igraph_vit_t tovit;
            IGRAPH_CHECK(igraph_vit_create(graph, to, &tovit));
            IGRAPH_FINALLY(igraph_vit_destroy, &tovit);
            for (j = 0, IGRAPH_VIT_RESET(tovit); j < nc; j++, IGRAPH_VIT_NEXT(tovit)) {
                long int v2 = IGRAPH_VIT_GET(tovit);
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
