/*
   IGraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_constructors.h"
#include "igraph_operators.h"

#include "igraph_conversion.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

/**
 * \function igraph_mycielskian
 * \brief Generate the Mycielskian of a graph with k iterations.
 *
 * The Mycielskian construction is a method used to generate a triangle-free
 * graph with an increased chromatic number. Given an input graph, the
 * Mycielskian transformation produces a new graph where the chromatic number
 * increases by at most one in each iteration.
 *
 * </para><para>
 * This transformation is useful in graph theory for constructing
 * graphs with large chromatic numbers while avoiding small cycles,
 * particularly triangles.
 *
 * \subsection mycielskian_structure Structure of Mycielskian Graph
 *
 * Let the \( n \) vertices of the given graph \( G \) be \( v_1, v_2, \dots, v_n \).
 * The Mycielski graph \( \mu(G) \) contains \( G \) itself as a subgraph,
 * together with \( n+1 \) additional vertices:
 * - A vertex \( u_i \) corresponding to each vertex \( v_i \) of \( G \).
 * - An extra vertex \( w \).
 *
 * The edges are added as follows:
 * - Each vertex \( u_i \) is connected to \( w \), forming a star \( K_{1,n} \).
 * - For each edge \( (v_i, v_j) \) in \( G \), two new edges are added:
 *   \( (u_i, v_j) \) and \( (v_i, u_j) \).
 *
 * Thus, if \( G \) has \( n \) vertices and \( m \) edges, the Mycielski graph \( \mu(G) \) has:
 * \[
 * 2n + 1 \text{ vertices, and } 3m + n \text{ edges.}
 * \]
 *
 * \subsection mycielskian_growth Growth with k Iterations
 * The k-th iterated Mycielskian has:
 * \f[
 * n_k = ( n + 1 ) \cdot 2^k - 1
 * \f]
 * vertices, where \( n \) is the number of vertices in the original graph.
 *
 * The edge count increases as:
 * \f[
 * m_k = \frac{(2m + 2n + 1) \cdot 3^k - n_{k+1}}{2}
 * \f]
 * where \( m \) is the original number of edges. The number of edges increases
 * from \( m \) to \( 3m + n \) in each iteration, as we add 2 edges for every
 * existing edge and 1 edge for every vertex.
 *
 * </para><para>
 * The Mycielskian is commonly used to demonstrate that graphs can have high
 * chromatic numbers without forming small cycles (i.e., triangle-free graphs).
 * The Grötzsch graph, for example, is obtained as \( M_4 \) of a triangle-free graph.
 *
 * \param graph Pointer to the input graph.
 * \param res Pointer to an uninitialized graph object where the Mycielskian
 *        of the input graph will be stored.
 * \param k Integer, the number of Mycielskian iterations (must be ≥ 0).
 * \return Error code.
 *
 * \sa \ref igraph_mycielski_graph().
 *
 * Time complexity: Exponential in k.
 */
igraph_error_t igraph_mycielskian(const igraph_t *graph, igraph_t *res, igraph_integer_t k) {
    if (k < 0) {
        IGRAPH_ERROR("The number of Mycielski iterations must not be negative.", IGRAPH_EINVAL);
    }

    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_integer_t ecount = igraph_ecount(graph);
    igraph_integer_t init_iter = 0; // to track the init iteration number
    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (vcount == 0) { // empty graph
        if (k <= 1) { // 0--> Null graph, 1--> single vertex
            IGRAPH_CHECK(igraph_empty(res, k, IGRAPH_UNDIRECTED));
            IGRAPH_FINALLY_CLEAN(1);
            igraph_vector_int_destroy(&edges);
            return IGRAPH_SUCCESS;
        }
        /* else */
        init_iter = 2; // start from 2nd iteration
        // Make a 2-vertex graph
        vcount = 2;
        ecount = 1;
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, 0));
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, 1));
    }

    if (vcount == 1) { // single vertex, assuming no self loop
        if (k == 0) {
            IGRAPH_CHECK(igraph_empty(res, 1, IGRAPH_UNDIRECTED));
            IGRAPH_FINALLY_CLEAN(1);
            igraph_vector_int_destroy(&edges);
            return IGRAPH_SUCCESS;
        }
        /* else */
        init_iter = 1; // start from 1st iteration
        // Make a 2-vertex graph
        vcount = 2;
        ecount = 1;
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, 0));
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, 1));
    }

    igraph_integer_t new_vcount = vcount;
    igraph_integer_t new_ecount = ecount;

    for (igraph_integer_t i = init_iter; i < k; i++) {
        // new edges = 3 * old edges + old vertices
        IGRAPH_SAFE_MULT(new_ecount, 3, &new_ecount);
        IGRAPH_SAFE_ADD(new_ecount, new_vcount, &new_ecount);

        // new vertices = 2 * old vertices + 1
        IGRAPH_SAFE_MULT(new_vcount, 2, &new_vcount);
        IGRAPH_SAFE_ADD(new_vcount, 1, &new_vcount);
    }

    // copy the edges from the original graph to the new vector
    // excluding for the case of vcount==0 and vcount==1
    // i.e, init_iter != 0, exclude it
    if (init_iter == 0) {
        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, false));
    }
    IGRAPH_CHECK(igraph_vector_int_resize(&edges, new_ecount * 2));

    igraph_integer_t edge_index = 2 * ecount;  // Current last edge index in edge vector
    igraph_integer_t offset = vcount;          // Tracks where new vertices start

    for (igraph_integer_t i = init_iter; i < k; i++) {
        igraph_integer_t prev_vcount = offset;  // Number of vertices before this step
        igraph_integer_t w = offset * 2;        // The new 'w' node index
        igraph_integer_t last_edge_index = edge_index;  // Mark where edges before this step end

        // For each edge before this step, add two new edges
        for (igraph_integer_t j = 0; j < last_edge_index; j += 2) {
            igraph_integer_t v1 = VECTOR(edges)[j];
            igraph_integer_t v2 = VECTOR(edges)[j + 1];

            VECTOR(edges)[edge_index++] = v1;
            VECTOR(edges)[edge_index++] = offset + v2;

            VECTOR(edges)[edge_index++] = v2;
            VECTOR(edges)[edge_index++] = offset + v1;
        }

        // Add edges connecting each `ui` to `w` (forming a star)
        for (igraph_integer_t j = prev_vcount; j < w; j++) {
            VECTOR(edges)[edge_index++] = j;
            VECTOR(edges)[edge_index++] = w;
        }

        // Update offset for next step
        offset = offset * 2 + 1;

        IGRAPH_ALLOW_INTERRUPTION();
    }

    // Add all edges in one go
    IGRAPH_CHECK(igraph_create(res, &edges, 0, igraph_is_directed(graph)));
    IGRAPH_FINALLY_CLEAN(1);
    igraph_vector_int_destroy(&edges);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_mycielski_graph
 * \brief Generate the Mycielski graph of order k.
 *
 * The Mycielski graph construction is used to create a triangle-free graph
 * with an increased chromatic number. Given an input graph, the Mycielski
 * transformation produces a new graph with one additional chromatic number.
 * This transformation is commonly used in graph theory to study chromatic
 * properties of graphs.
 *
 * </para><para>
 * The Mycielski graph \( M_k \) of order \( k \) is a triangle-free graph
 * with chromatic number \( k \) and the smallest possible number of vertices.
 *
 * The number of edges also increases with each step, making it useful for
 * demonstrating properties of chromatic graphs.
 *
 * The first few Mycielski graphs are:
 * - \( M_0 \): Null graph
 * - \( M_1 \): Single vertex
 * - \( M_2 \): Path graph with 2 vertices
 * - \( M_3 \): Cycle graph with 5 vertices
 * - \( M_4 \): Grötzsch graph (triangle-free graph with chromatic number 4)
 *
 * </para><para>
 * The Mycielski graph has several important applications, particularly in
 * demonstrating that there exist graphs with high chromatic numbers but no
 * short cycles (i.e., triangle-free). For example, the Grötzsch graph,
 * which is \( M_4 \), is a triangle-free graph with chromatic number 4.
 *
 * \param graph Pointer to an uninitialized graph object. The generated
 *        Mycielski graph will be stored here.
 * \param k Integer, the order of the Mycielski graph (must be \( \geq 0 \)).
 * \return Error code.
 *
 * \sa \ref igraph_mycielskian().
 *
 * Time complexity: Exponential in \( k \).
 */
igraph_error_t igraph_mycielski_graph(igraph_t *graph, igraph_integer_t k) {
    igraph_t g;

    if (k < 0) {
        IGRAPH_ERROR("The Mycielski graph order must not be negative.", IGRAPH_EINVAL);
    }
    if (k <= 1) { // 0--> Null graph, 1--> single vertex
        IGRAPH_CHECK(igraph_empty(graph, k, IGRAPH_UNDIRECTED));
        return IGRAPH_SUCCESS;
    }

    // create a path 0---1
    IGRAPH_CHECK(igraph_ring(&g, 2, IGRAPH_UNDIRECTED, false, false));
    IGRAPH_FINALLY(igraph_destroy, &g);

    IGRAPH_CHECK(igraph_mycielskian(&g, graph, k - 2));

    igraph_destroy(&g);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
