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

#include "igraph_operators.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"

#include "math/safe_intop.h"

static igraph_error_t cartesian_product(igraph_t *res,
                                        const igraph_t *g1,
                                        const igraph_t *g2) {

    const igraph_bool_t directed = igraph_is_directed(g1);

    if (igraph_is_directed(g2) != directed) {
        IGRAPH_ERROR("Cartesian product between a directed and an undirected graph is invalid.",
                     IGRAPH_EINVAL);
    }

    const igraph_integer_t vcount1 = igraph_vcount(g1);
    const igraph_integer_t vcount2 = igraph_vcount(g2);
    const igraph_integer_t ecount1 = igraph_ecount(g1);
    const igraph_integer_t ecount2 = igraph_ecount(g2);
    igraph_integer_t vcount;
    igraph_integer_t ecount, ecount_double;
    igraph_vector_int_t edges;

    // New vertex count = vcount1 * vcount2
    IGRAPH_SAFE_MULT(vcount1, vcount2, &vcount);

    {
        // New edge count = vcount1*ecount2 + vcount2*ecount1
        igraph_integer_t temp;
        IGRAPH_SAFE_MULT(vcount1, ecount2, &ecount);
        IGRAPH_SAFE_MULT(vcount2, ecount1, &temp);
        IGRAPH_SAFE_ADD(ecount, temp, &ecount);
    }

    IGRAPH_SAFE_MULT(ecount, 2, &ecount_double);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ecount_double);

    // Vertex ((i, j)) with i from g1, and j from g2
    //   will have new vertex id: i * vcount2 + j

    igraph_integer_t edge_index = 0;

    // Edges from g1
    for (igraph_integer_t i = 0; i < ecount1; ++i) {
        igraph_integer_t from = IGRAPH_FROM(g1, i);
        igraph_integer_t to = IGRAPH_TO(g1, i);

        // For all edges (from, to) in g1, add edge from ((from, j)) to ((to, j))
        //    for all vertex j in g2
        for (igraph_integer_t j = 0; j < vcount2; ++j) {
            VECTOR(edges)[edge_index++] = from * vcount2 + j; // ((from, j))
            VECTOR(edges)[edge_index++] = to * vcount2 + j; // ((to, j))
        }
    }

    // Edges from g2
    for (igraph_integer_t i = 0; i < ecount2; ++i) {
        igraph_integer_t from = IGRAPH_FROM(g2, i);
        igraph_integer_t to = IGRAPH_TO(g2, i);

        // For all edges (from, to) in g2, add edge from (j, from) to (j, to)
        //    for all vertex j in g1
        for (igraph_integer_t j = 0; j < vcount1; ++j) {
            VECTOR(edges)[edge_index++] = j * vcount2 + from; // ((j, from))
            VECTOR(edges)[edge_index++] = j * vcount2 + to; // ((j, to))
        }
    }

    IGRAPH_CHECK(igraph_create(res, &edges, vcount, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t lexicographic_product(igraph_t *res,
                                            const igraph_t *g1,
                                            const igraph_t *g2) {

    const igraph_bool_t directed = igraph_is_directed(g1);

    if (igraph_is_directed(g2) != directed) {
        IGRAPH_ERROR("Lexicographic product between a directed and an undirected graph is invalid.",
                     IGRAPH_EINVAL);
    }

    const igraph_integer_t vcount1 = igraph_vcount(g1);
    const igraph_integer_t vcount2 = igraph_vcount(g2);
    const igraph_integer_t ecount1 = igraph_ecount(g1);
    const igraph_integer_t ecount2 = igraph_ecount(g2);
    igraph_integer_t vcount;
    igraph_integer_t ecount, ecount_double;
    igraph_vector_int_t edges;

    // New vertex count = vcount1 * vcount2
    IGRAPH_SAFE_MULT(vcount1, vcount2, &vcount);

    {
        // New edge count = vcount1*ecount2 + (vcount2^2)*ecount1
        igraph_integer_t temp;
        IGRAPH_SAFE_MULT(vcount1, ecount2, &ecount);
        IGRAPH_SAFE_MULT(vcount2, vcount2, &temp);
        IGRAPH_SAFE_MULT(temp, ecount1, &temp);

        IGRAPH_SAFE_ADD(ecount, temp, &ecount);
    }

    IGRAPH_SAFE_MULT(ecount, 2, &ecount_double);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ecount_double);

    // Vertex ((i, j)) with i from g1, and j from g2
    //   will have new vertex id: i * vcount2 + j
    igraph_integer_t edge_index = 0;

    // edges of form u1=u2 and v1~v2
    for (igraph_integer_t i = 0; i < ecount2; ++i) {
        igraph_integer_t from = IGRAPH_FROM(g2, i);
        igraph_integer_t to = IGRAPH_TO(g2, i);

        // For all edges (from, to) in g2, add edge from (j, from) to (j, to)
        //    for all vertex j in g1
        for (igraph_integer_t j = 0; j < vcount1; ++j) {
            VECTOR(edges)[edge_index++] = j * vcount2 + from; // ((j, from))
            VECTOR(edges)[edge_index++] = j * vcount2 + to; // ((j, to))
        }
    }

    // edges of form u1~u2
    for (igraph_integer_t i = 0; i < ecount1; ++i) {
        igraph_integer_t from1 = IGRAPH_FROM(g1, i);
        igraph_integer_t to1 = IGRAPH_TO(g1, i);

        // each vertex pair irrespective of their connectivity
        for (igraph_integer_t from2 = 0; from2 < vcount2; ++from2) {
            for (igraph_integer_t to2 = 0; to2 < vcount2; ++to2) {
                // ((from1, from2)) to ((to1, to2))
                VECTOR(edges)[edge_index++] = from1 * vcount2 + from2; // ((from1, from2))
                VECTOR(edges)[edge_index++] = to1 * vcount2 + to2; // ((to1, to2))
            }
        }
    }

    IGRAPH_CHECK(igraph_create(res, &edges, vcount, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t strong_product(igraph_t *res,
                                     const igraph_t *g1,
                                     const igraph_t *g2) {

    const igraph_bool_t directed = igraph_is_directed(g1);

    if (igraph_is_directed(g2) != directed) {
        IGRAPH_ERROR("Strong product between a directed and an undirected graph is invalid.",
                     IGRAPH_EINVAL);
    }

    const igraph_integer_t vcount1 = igraph_vcount(g1);
    const igraph_integer_t vcount2 = igraph_vcount(g2);
    const igraph_integer_t ecount1 = igraph_ecount(g1);
    const igraph_integer_t ecount2 = igraph_ecount(g2);
    igraph_integer_t vcount;
    igraph_integer_t ecount, ecount_double;
    igraph_vector_int_t edges;

    // New vertex count = vcount1 * vcount2
    IGRAPH_SAFE_MULT(vcount1, vcount2, &vcount);

    {
        // New edge count = vcount1*ecount2 + vcount2*ecount1 + 2*e1*e2 for undirected graph
        //                = vcount1*ecount2 + vcount2*ecount1 + e1*e2 for directed graph
        igraph_integer_t temp;
        IGRAPH_SAFE_MULT(vcount1, ecount2, &ecount);
        IGRAPH_SAFE_MULT(vcount2, ecount1, &temp);
        IGRAPH_SAFE_ADD(ecount, temp, &ecount);

        IGRAPH_SAFE_MULT(ecount1, ecount2, &temp);
        if (!directed) {
            IGRAPH_SAFE_MULT(temp, 2, &temp);
        }
        IGRAPH_SAFE_ADD(ecount, temp, &ecount);
    }

    IGRAPH_SAFE_MULT(ecount, 2, &ecount_double);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ecount_double);

    igraph_integer_t edge_index = 0;
    // Strong graph product contains all the edges from both cartesian and tensor product

    // Edges of type cartesian products: v1=v2 and u1~u2; v1~v2 and u1=u2
    // v1=v2 and u1~u2
    for (igraph_integer_t i = 0; i < ecount1; ++i) {
        igraph_integer_t from = IGRAPH_FROM(g1, i);
        igraph_integer_t to = IGRAPH_TO(g1, i);

        // For all edges (from, to) in g1, add edge from ((from, j)) to ((to, j))
        //    for all vertex j in g2
        for (igraph_integer_t j = 0; j < vcount2; ++j) {
            // SAFE MULT and SAFE ADD not needed as < vcount
            VECTOR(edges)[edge_index++] = from * vcount2 + j; // ((from, j))
            VECTOR(edges)[edge_index++] = to * vcount2 + j; // ((to, j))
        }
    }

    // v1~v2 and u1=u2
    for (igraph_integer_t i = 0; i < ecount2; ++i) {
        igraph_integer_t from = IGRAPH_FROM(g2, i);
        igraph_integer_t to = IGRAPH_TO(g2, i);

        // For all edges (from, to) in g2, add edge from (j, from) to (j, to)
        //    for all vertex j in g1
        for (igraph_integer_t j = 0; j < vcount1; ++j) {
            VECTOR(edges)[edge_index++] = j * vcount2 + from; // ((j, from))
            VECTOR(edges)[edge_index++] = j * vcount2 + to; // ((j, to))
        }
    }

    // Edges of type tensor product
    // u1 ~ u2 and v1 ~ v2
    for (igraph_integer_t i = 0; i < ecount1; ++i) {
        igraph_integer_t from1 = IGRAPH_FROM(g1, i);
        igraph_integer_t to1 = IGRAPH_TO(g1, i);

        for (igraph_integer_t j = 0; j < ecount2; ++j) {
            igraph_integer_t from2 = IGRAPH_FROM(g2, j);
            igraph_integer_t to2 = IGRAPH_TO(g2, j);

            // Create edge between ((from1, from2)) to ((to1, to2))
            VECTOR(edges)[edge_index++] = from1 * vcount2 + from2; // ((from1, from2))
            VECTOR(edges)[edge_index++] = to1 * vcount2 + to2; // ((to1, to2))

            // In directed graphs, no edge is added because (from2, to2) are not adjacent
            // respecting direction.
            // For undirected graphs, add cross edges between ((from1, to2)) and ((to1, from2)).
            if (!directed) {
                VECTOR(edges)[edge_index++] = from1 * vcount2 + to2; // ((from1, to2))
                VECTOR(edges)[edge_index++] = to1 * vcount2 + from2; // ((to1, from2))
            }
        }
    }
    IGRAPH_CHECK(igraph_create(res, &edges, vcount, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t tensor_product(igraph_t *res,
                                     const igraph_t *g1,
                                     const igraph_t *g2) {

    const igraph_bool_t directed = igraph_is_directed(g1);

    if (igraph_is_directed(g2) != directed) {
        IGRAPH_ERROR("Tensor product between a directed and an undirected graph is invalid.",
                     IGRAPH_EINVAL);
    }

    const igraph_integer_t vcount1 = igraph_vcount(g1);
    const igraph_integer_t vcount2 = igraph_vcount(g2);
    const igraph_integer_t ecount1 = igraph_ecount(g1);
    const igraph_integer_t ecount2 = igraph_ecount(g2);
    igraph_integer_t vcount;
    igraph_integer_t ecount, ecount_double;
    igraph_vector_int_t edges;

    IGRAPH_SAFE_MULT(vcount1, vcount2, &vcount);

    // New edge count = 2*ecount1*ecount2 if undirected else ecount1*ecount2
    IGRAPH_SAFE_MULT(ecount1, ecount2, &ecount);
    if (!directed) {
        IGRAPH_SAFE_MULT(ecount, 2, &ecount);
    }
    IGRAPH_SAFE_MULT(ecount, 2, &ecount_double);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ecount_double);

    // Vertex ((i, j)) with i from g1, and j from g2
    //   will have new vertex id: i * vcount2 + j

    igraph_integer_t edge_index = 0;

    for (igraph_integer_t i = 0; i < ecount1; ++i) {
        igraph_integer_t from1 = IGRAPH_FROM(g1, i);
        igraph_integer_t to1 = IGRAPH_TO(g1, i);

        for (igraph_integer_t j = 0; j < ecount2; ++j) {
            igraph_integer_t from2 = IGRAPH_FROM(g2, j);
            igraph_integer_t to2 = IGRAPH_TO(g2, j);

            // Create edge between ((from1, from2)) to ((to1, to2))
            VECTOR(edges)[edge_index++] = from1 * vcount2 + from2; // ((from1, from2))
            VECTOR(edges)[edge_index++] = to1 * vcount2 + to2; // ((to1, to2))

            // In directed graphs, no edge is added because (from2, to2) are not adjacent
            // respecting direction.
            // For undirected graphs, add cross edges between ((from1, to2)) and ((to1, from2)).
            if (!directed) {
                VECTOR(edges)[edge_index++] = from1 * vcount2 + to2; // ((from1, to2))
                VECTOR(edges)[edge_index++] = to1 * vcount2 + from2; // ((to1, from2))
            }
        }
    }

    IGRAPH_CHECK(igraph_create(res, &edges, vcount, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_product
 * \brief The graph product of two graphs, according to the chosen product type.
 *
 * \experimental
 *
 * This function computes the product of two graphs using the graph
 * product concept selected by the \p type parameter. The two graphs must be of
 * the same type, either directed or undirected. If a product of an undirected
 * and a directed graph is required, convert one of them to the appropriate type
 * using \ref igraph_to_directed() or \ref igraph_to_undirected().
 *
 * </para><para>
 * Each vertex of the product graph corresponds to a pair <code>(u, v)</code>,
 * where \c u is a vertex from the first graph and \c v is a vertex from the
 * second graph. Thus the number of vertices in the product graph is
 * <code>|V1| |V2|</code>,
 * where <code>|V1|</code> and <code>|V2|</code> are the sizes of the vertex
 * set of the operands. The pair <code>(u, v)</code> is mapped to a unique vertex
 * index in the product graph using <code>index = u |V2| + v</code>.
 *
 * </para><para>
 * All implemented graph products are associative, but not all are commutative.
 *
 * </para><para>
 * The supported graph product types are detailed below. The notation
 * <code>u ~ v</code>
 * is used to indicate that vertices \c u and \c v are adjacent, i.e. there is
 * a connection from \c u to \c v.
 * \clist
 *     \cli IGRAPH_PRODUCT_CARTESIAN
 *     Computes the cartesian product of two graphs. In the product graph,
 *     there is a connection from <code>(u1, v1)</code> to <code>(u2, v2)</code>
 *     if and only if
 *     <code>u1 = u2</code> and <code>v1 ~ v2</code> or
 *     <code>u1 ~ u2</code> and <code>v1 = v2</code>.
 *     Thus, the number of edges in the product graph is
 *     <code>|V1| |E2| + |V2| |E1|</code>.
 *
 *     </para><para>
 *     Time complexity: O(|V1| |V2| + |V1| |E2| + |V2| |E1|)
 *     where |V1| and |V2| are the number of vertices, and
 *     |E1| and |E2| are the number of edges of the operands.
 *
 *     \cli IGRAPH_PRODUCT_LEXICOGRAPHIC
 *     Computes the lexicographic product of two graphs. In the product graph,
 *     there is a connection from <code>(u1, v1)</code> to <code>(u2, v2)</code>
 *     if and only if
 *     <code>u1 = u2</code> and <code>v1 ~ v2</code> or
 *     <code>u1 ~ u2</code>.
 *     Thus, the number of edges in the product graph is
 *     <code>|V1| |E2| + |V2|^2 |E1|</code>. Unlike most other graph products,
 *     the lexicographic product is not commutative.
 *
 *     </para><para>
 *     Time complexity: O(|V1| |V2| + |V1| |E2| + |V2|^2 |E1|)
 *     where |V1| and |V2| are the number of vertices, and
 *     |E1| and |E2| are the number of edges of the operands.
 *
 *     \cli IGRAPH_PRODUCT_STRONG
 *     Computes the strong product (also known as normal product) of two graphs.
 *     In the product graph, there is a connection from <code>(u1, v1)</code> to
 *     <code>(u2, v2)</code>
 *     if and only if
 *     <code>u1 = u2</code> and <code>v1 ~ v2</code> or
 *     <code>u1 ~ u2</code> and <code>v1 = v2</code> or
 *     <code>u1 ~ u2</code> and <code>v1 ~ v2</code>.
 *     Thus, the number of edges in the product graph is
 *     <code>|V1| |E2| + |V2| |E1| + |E1| |E2|</code> in the directed case and
 *     <code>|V1| |E2| + |V2| |E1| + 2 |E1| |E2|</code> in the undirected case.
 *
 *     </para><para>
 *     Time complexity: O(|V1| |V2| + |V1| |E2| + |V2| |E1| + |E1| |E2|)
 *     where |V1| and |V2| are the number of vertices, and
 *     |E1| and |E2| are the number of edges of the operands.
 *
 *     \cli IGRAPH_PRODUCT_TENSOR
 *     Computes the tensor product (also known as categorial product) of two graphs.
 *     In the product graph, there is a connection from <code>(u1, v1)</code> to
 *     <code>(u2, v2)</code>
 *     if and only if
 *     <code>u1 ~ u2</code> and <code>v1 ~ v2</code>.
 *     Thus, the number of edges in the product is
 *     <code>|E1| |E2|</code> in the directed case and
 *     <code>2 |E1| |E2|</code> in the undirected case.
 *
 *     </para><para>
 *     Time Complexity: O(|V1| |V2| + |E1| |E2|)
 *     where |V1| and |V2| are the number of vertices, and
 *     |E1| and |E2| are the number of edges of the operands.
 * \endclist
 *
 * </para><para>
 * Reference:
 *
  * </para><para>
 * Richar Hammack, Wilfried Imrich, and Sandi Klavžar (2011).
 * Handbook of Product Graphs (2nd ed.). CRC Press.
 * https://doi.org/10.1201/b10959
 *
 * \param res Pointer to an uninitialized graph object. The product graph will
 *   be stored here.
 * \param g1 The first operand graph.
 * \param g2 The second operand graph. It must have the same directedness as \p g1.
 * \param type The type of graph product to compute.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL if the specified \p type is unsupported or the input
 *         graphs \p g1 and \p g2 are incompatible for the requested product.
 */

igraph_error_t igraph_product(igraph_t *res,
                              const igraph_t *g1,
                              const igraph_t *g2,
                              igraph_product_t type) {
    switch (type) {
    case IGRAPH_PRODUCT_CARTESIAN:
        return cartesian_product(res, g1, g2);

    case IGRAPH_PRODUCT_LEXICOGRAPHIC:
        return lexicographic_product(res, g1, g2);

    case IGRAPH_PRODUCT_STRONG:
        return strong_product(res, g1, g2);

    case IGRAPH_PRODUCT_TENSOR:
        return tensor_product(res, g1, g2);

    default:
        IGRAPH_ERROR("Unknown graph product type.", IGRAPH_EINVAL);
    }
}
