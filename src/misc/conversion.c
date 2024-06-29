/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_conversion.h"

#include "igraph_iterators.h"
#include "igraph_interface.h"
#include "igraph_attributes.h"
#include "igraph_constructors.h"
#include "igraph_structural.h"
#include "igraph_sparsemat.h"
#include "igraph_random.h"

#include "core/fixed_vectorlist.h"
#include "graph/attributes.h"
#include "math/safe_intop.h"

#define WEIGHT_OF(eid) (weights ? VECTOR(*weights)[eid] : 1)

/**
 * \ingroup conversion
 * \function igraph_get_adjacency
 * \brief The adjacency matrix of a graph.
 *
 * </para><para>
 * The result is an adjacency matrix. Entry i, j of the matrix
 * contains the number of edges connecting vertex i to vertex j in the unweighted
 * case, or the total weight of edges connecting vertex i to vertex j in the
 * weighted case.
 *
 * \param graph Pointer to the graph to convert
 * \param res Pointer to an initialized matrix object, it will be
 *        resized if needed.
 * \param type Constant specifying the type of the adjacency matrix to
 *        create for undirected graphs. It is ignored for directed
 *        graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_GET_ADJACENCY_UPPER
 *          the upper right triangle of the matrix is used.
 *        \cli IGRAPH_GET_ADJACENCY_LOWER
 *          the lower left triangle of the matrix is used.
 *        \cli IGRAPH_GET_ADJACENCY_BOTH
 *          the whole matrix is used, a symmetric matrix is returned
 *          if the graph is undirected.
 *        \endclist
 * \param weights An optional vector containing the weight of each edge
 *        in the graph. Supply a null pointer here to make all edges have
 *        the same weight of 1.
 * \param loops Constant specifying how loop edges should be handled.
 *        Possible values:
 *        \clist
 *        \cli IGRAPH_NO_LOOPS
 *          loop edges are ignored and the diagonal of the matrix will contain
 *          zeros only
 *        \cli IGRAPH_LOOPS_ONCE
 *          loop edges are counted once, i.e. a vertex with a single unweighted
 *          loop edge will have 1 in the corresponding diagonal entry
 *        \cli IGRAPH_LOOPS_TWICE
 *          loop edges are counted twice in \em undirected graphs, i.e. a vertex
 *          with a single unweighted loop edge in an undirected graph will have
 *          2 in the corresponding diagonal entry. Loop edges in directed graphs
 *          are still counted as 1. Essentially, this means that the function is
 *          counting the incident edge \em stems , which makes more sense when
 *          using the adjacency matrix in linear algebra.
 *        \endclist
 * \return Error code:
 *        \c IGRAPH_EINVAL invalid type argument.
 *
 * \sa \ref igraph_get_adjacency_sparse() if you want a sparse matrix representation
 *
 * Time complexity: O(|V||V|), |V| is the number of vertices in the graph.
 */

igraph_error_t igraph_get_adjacency(
    const igraph_t *graph, igraph_matrix_t *res, igraph_get_adjacency_t type,
    const igraph_vector_t *weights, igraph_loops_t loops
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t i, from, to;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);

    if (directed) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (from != to || loops != IGRAPH_NO_LOOPS) {
                MATRIX(*res, from, to) += WEIGHT_OF(i);
            }
        }
    } else if (type == IGRAPH_GET_ADJACENCY_UPPER) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (to < from) {
                MATRIX(*res, to, from) += WEIGHT_OF(i);
            } else {
                MATRIX(*res, from, to) += WEIGHT_OF(i);
            }
            if (to == from && loops == IGRAPH_LOOPS_TWICE) {
                MATRIX(*res, to, to) += WEIGHT_OF(i);
            }
        }
    } else if (type == IGRAPH_GET_ADJACENCY_LOWER) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (to < from) {
                MATRIX(*res, from, to) += WEIGHT_OF(i);
            } else {
                MATRIX(*res, to, from) += WEIGHT_OF(i);
            }
            if (to == from && loops == IGRAPH_LOOPS_TWICE) {
                MATRIX(*res, to, to) += WEIGHT_OF(i);
            }
        }
    } else if (type == IGRAPH_GET_ADJACENCY_BOTH) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            MATRIX(*res, from, to) += WEIGHT_OF(i);
            if (from != to || loops == IGRAPH_LOOPS_TWICE) {
                MATRIX(*res, to, from) += WEIGHT_OF(i);
            }
        }
    } else {
        IGRAPH_ERROR("Invalid type argument", IGRAPH_EINVAL);
    }

    /* Erase the diagonal if we don't need loop edges */
    if (loops == IGRAPH_NO_LOOPS) {
        for (i = 0; i < no_of_nodes; i++) {
            MATRIX(*res, i, i) = 0;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_adjacency_sparse
 * \brief Returns the adjacency matrix of a graph in a sparse matrix format.
 *
 * \param graph The input graph.
 * \param res Pointer to an \em initialized sparse matrix. The result
 *    will be stored here. The matrix will be resized as needed.
 * \param type Constant specifying the type of the adjacency matrix to
 *        create for undirected graphs. It is ignored for directed
 *        graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_GET_ADJACENCY_UPPER
 *          the upper right triangle of the matrix is used.
 *        \cli IGRAPH_GET_ADJACENCY_LOWER
 *          the lower left triangle of the matrix is used.
 *        \cli IGRAPH_GET_ADJACENCY_BOTH
 *          the whole matrix is used, a symmetric matrix is returned
 *          if the graph is undirected.
 *        \endclist
 * \return Error code:
 *        \c IGRAPH_EINVAL invalid type argument.
 *
 * \sa \ref igraph_get_adjacency(), the dense version of this function.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_get_adjacency_sparse(
    const igraph_t *graph, igraph_sparsemat_t *res, igraph_get_adjacency_t type,
    const igraph_vector_t *weights, igraph_loops_t loops
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t nzmax = directed ? no_of_edges : no_of_edges * 2;
    igraph_integer_t i, from, to;

    IGRAPH_CHECK(igraph_sparsemat_resize(res, no_of_nodes, no_of_nodes, nzmax));

    if (directed) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (from != to || loops != IGRAPH_NO_LOOPS) {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, from, to, WEIGHT_OF(i)));
            }
        }
    } else if (type == IGRAPH_GET_ADJACENCY_UPPER) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (to < from) {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, to, from, WEIGHT_OF(i)));
            } else if (to == from) {
                switch (loops) {
                    case IGRAPH_LOOPS_ONCE:
                        IGRAPH_CHECK(igraph_sparsemat_entry(res, to, to, WEIGHT_OF(i)));
                        break;
                    case IGRAPH_LOOPS_TWICE:
                        IGRAPH_CHECK(igraph_sparsemat_entry(res, to, to, 2 * WEIGHT_OF(i)));
                        break;
                    case IGRAPH_NO_LOOPS:
                    default:
                        break;
                }
            } else {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, from, to, WEIGHT_OF(i)));
            }
        }
    } else if (type == IGRAPH_GET_ADJACENCY_LOWER) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (to < from) {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, from, to, WEIGHT_OF(i)));
            } else if (to == from) {
                switch (loops) {
                    case IGRAPH_LOOPS_ONCE:
                        IGRAPH_CHECK(igraph_sparsemat_entry(res, to, to, WEIGHT_OF(i)));
                        break;
                    case IGRAPH_LOOPS_TWICE:
                        IGRAPH_CHECK(igraph_sparsemat_entry(res, to, to, 2 * WEIGHT_OF(i)));
                        break;
                    case IGRAPH_NO_LOOPS:
                    default:
                        break;
                }
            } else {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, to, from, WEIGHT_OF(i)));
            }
        }
    } else if (type == IGRAPH_GET_ADJACENCY_BOTH) {
        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            if (to == from) {
                switch (loops) {
                    case IGRAPH_LOOPS_ONCE:
                        IGRAPH_CHECK(igraph_sparsemat_entry(res, to, to, WEIGHT_OF(i)));
                        break;
                    case IGRAPH_LOOPS_TWICE:
                        IGRAPH_CHECK(igraph_sparsemat_entry(res, to, to, 2 * WEIGHT_OF(i)));
                        break;
                    case IGRAPH_NO_LOOPS:
                    default:
                        break;
                }
            } else {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, from, to, WEIGHT_OF(i)));
                IGRAPH_CHECK(igraph_sparsemat_entry(res, to, from, WEIGHT_OF(i)));
            }
        }
    } else {
        IGRAPH_ERROR("Invalid type argument", IGRAPH_EINVAL);
    }

    return IGRAPH_SUCCESS;
}

#undef WEIGHT_OF

/**
 * \function igraph_get_sparsemat
 * \brief Converts an igraph graph to a sparse matrix (deprecated).
 *
 * If the graph is undirected, then a symmetric matrix is created.
 *
 * </para><para>
 * This function is deprecated in favour of \ref igraph_get_adjacency_sparse(),
 * but does not work in an identical way. This function takes an \em uninitialized
 * \c igraph_sparsemat_t while \ref igraph_get_adjacency_sparse() takes
 * an already initialized one.
 *
 * \param graph The input graph.
 * \param res Pointer to an \em uninitialized sparse matrix. The result
 *    will be stored here.
 * \return Error code.
 *
 * \deprecated-by igraph_get_adjacency_sparse 0.10.0
 */

igraph_error_t igraph_get_sparsemat(const igraph_t *graph, igraph_sparsemat_t *res) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t nzmax = igraph_is_directed(graph) ? no_of_edges : 2*no_of_edges;
    IGRAPH_CHECK(igraph_sparsemat_init(res, no_of_nodes, no_of_nodes, nzmax));
    return igraph_get_adjacency_sparse(graph, res, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_ONCE);
}

/**
 * \ingroup conversion
 * \function igraph_get_edgelist
 * \brief The list of edges in a graph.
 *
 * The order of the edges is given by the edge IDs.
 *
 * \param graph Pointer to the graph object
 * \param res Pointer to an initialized vector object, it will be
 *        resized.
 * \param bycol Logical, if true, the edges will be returned
 *        columnwise, e.g. the first edge is
 *        <code>res[0]->res[|E|]</code>, the second is
 *        <code>res[1]->res[|E|+1]</code>, etc.
 * \return Error code.
 *
 * \sa \ref igraph_edges() to return the result only for some edge IDs.
 *
 * Time complexity: O(|E|), the number of edges in the graph.
 */

igraph_error_t igraph_get_edgelist(const igraph_t *graph, igraph_vector_int_t *res, igraph_bool_t bycol) {

    igraph_eit_t edgeit;
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t vptr = 0;
    igraph_integer_t from, to;

    IGRAPH_CHECK(igraph_vector_int_resize(res, no_of_edges * 2));
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                                   &edgeit));
    IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);

    if (bycol) {
        while (!IGRAPH_EIT_END(edgeit)) {
            igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &from, &to);
            VECTOR(*res)[vptr] = from;
            VECTOR(*res)[vptr + no_of_edges] = to;
            vptr++;
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else {
        while (!IGRAPH_EIT_END(edgeit)) {
            igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &from, &to);
            VECTOR(*res)[vptr++] = from;
            VECTOR(*res)[vptr++] = to;
            IGRAPH_EIT_NEXT(edgeit);
        }
    }

    igraph_eit_destroy(&edgeit);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_to_directed
 * \brief Convert an undirected graph to a directed one.
 *
 * If the supplied graph is directed, this function does nothing.
 *
 * \param graph The graph object to convert.
 * \param mode Constant, specifies the details of how exactly the
 *        conversion is done. Possible values:
 *        \clist
 *        \cli IGRAPH_TO_DIRECTED_ARBITRARY
 *        The number of edges in the
 *        graph stays the same, an arbitrarily directed edge is
 *        created for each undirected edge.
 *        \cli IGRAPH_TO_DIRECTED_MUTUAL
 *        Two directed edges are
 *        created for each undirected edge, one in each direction.
 *        \cli IGRAPH_TO_DIRECTED_RANDOM
 *        Each undirected edge is converted to a randomly oriented
 *        directed one.
 *        \cli IGRAPH_TO_DIRECTED_ACYCLIC
 *        Each undirected edge is converted to a directed edge oriented
 *        from a lower index vertex to a higher index one. If no self-loops
 *        were present, then the result is a directed acyclic graph.
 *        \endclist
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */

igraph_error_t igraph_to_directed(igraph_t *graph,
                       igraph_to_directed_t mode) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (igraph_is_directed(graph)) {
        return IGRAPH_SUCCESS;
    }

    switch (mode) {
    case IGRAPH_TO_DIRECTED_ARBITRARY:
    case IGRAPH_TO_DIRECTED_RANDOM:
    case IGRAPH_TO_DIRECTED_ACYCLIC:
      {
        igraph_t newgraph;
        igraph_vector_int_t edges;
        igraph_integer_t size = no_of_edges * 2;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, size);
        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));

        if (mode == IGRAPH_TO_DIRECTED_RANDOM) {
            RNG_BEGIN();

            for (igraph_integer_t i=0; i < no_of_edges; ++i) {
                if (RNG_INTEGER(0,1)) {
                    igraph_integer_t temp = VECTOR(edges)[2*i];
                    VECTOR(edges)[2*i] = VECTOR(edges)[2*i+1];
                    VECTOR(edges)[2*i+1] = temp;
                }
            }

            RNG_END();
        } else if (mode == IGRAPH_TO_DIRECTED_ACYCLIC) {
            /* Currently, the endpoints of undirected edges are ordered in the
               internal graph datastructure, i.e. it is always true that from < to.
               However, it is not guaranteed that this will not be changed in
               the future, and this ordering should not be relied on outside of
               the implementation of the minimal API in type_indexededgelist.c.

               Therefore, we order the edge endpoints anyway in the following loop: */
            for (igraph_integer_t i=0; i < no_of_edges; ++i) {
                if (VECTOR(edges)[2*i] > VECTOR(edges)[2*i+1]) {
                    igraph_integer_t temp = VECTOR(edges)[2*i];
                    VECTOR(edges)[2*i] = VECTOR(edges)[2*i+1];
                    VECTOR(edges)[2*i+1] = temp;
                }
            }
        }

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   no_of_nodes,
                                   IGRAPH_DIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, true, true, true);
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(2);

        igraph_destroy(graph);
        *graph = newgraph;

        break;
      }
    case IGRAPH_TO_DIRECTED_MUTUAL:
      {
        igraph_t newgraph;
        igraph_vector_int_t edges;
        igraph_vector_int_t index;
        igraph_integer_t size;

        IGRAPH_SAFE_MULT(no_of_edges, 4, &size);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, size));
        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
        IGRAPH_CHECK(igraph_vector_int_resize(&edges, size));
        IGRAPH_VECTOR_INT_INIT_FINALLY(&index, no_of_edges * 2);
        for (igraph_integer_t i = 0; i < no_of_edges; i++) {
            VECTOR(edges)[no_of_edges * 2 + i * 2]  = VECTOR(edges)[i * 2 + 1];
            VECTOR(edges)[no_of_edges * 2 + i * 2 + 1] = VECTOR(edges)[i * 2];
            VECTOR(index)[i] = VECTOR(index)[no_of_edges + i] = i;
        }

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   no_of_nodes,
                                   IGRAPH_DIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, true, true, /*edges=*/false);
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, &newgraph, &index));

        igraph_vector_int_destroy(&index);
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(3);

        igraph_destroy(graph);
        *graph = newgraph;

        break;
      }
    default:
        IGRAPH_ERROR("Cannot direct graph, invalid mode.", IGRAPH_EINVAL);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_to_undirected
 * \brief Convert a directed graph to an undirected one.
 *
 * </para><para>
 * If the supplied graph is undirected, this function does nothing.
 *
 * \param graph The graph object to convert.
 * \param mode Constant, specifies the details of how exactly the
 *        conversion is done. Possible values: \c
 *        IGRAPH_TO_UNDIRECTED_EACH: the number of edges remains
 *        constant, an undirected edge is created for each directed
 *        one, this version might create graphs with multiple edges;
 *        \c IGRAPH_TO_UNDIRECTED_COLLAPSE: one undirected edge will
 *        be created for each pair of vertices that are connected
 *        with at least one directed edge, no multiple edges will be
 *        created. \c IGRAPH_TO_UNDIRECTED_MUTUAL creates an undirected
 *        edge for each pair of mutual edges in the directed graph.
 *        Non-mutual edges are lost; loop edges are kept unconditionally.
 *        This mode might create multiple edges.
 * \param edge_comb What to do with the edge attributes. See the igraph
 *        manual section about attributes for details. \c NULL means that
 *        the edge attributes are lost during the conversion, \em except
 *        when \c mode is \c IGRAPH_TO_UNDIRECTED_EACH, in which case the
 *        edge attributes are kept intact.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 *
 * \example examples/simple/igraph_to_undirected.c
 */

igraph_error_t igraph_to_undirected(igraph_t *graph,
                         igraph_to_undirected_t mode,
                         const igraph_attribute_combination_t *edge_comb) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t edges;
    igraph_t newgraph;
    igraph_bool_t attr = edge_comb && igraph_has_attribute_table();

    if (mode != IGRAPH_TO_UNDIRECTED_EACH &&
        mode != IGRAPH_TO_UNDIRECTED_COLLAPSE &&
        mode != IGRAPH_TO_UNDIRECTED_MUTUAL) {
        IGRAPH_ERROR("Cannot undirect graph, invalid mode", IGRAPH_EINVAL);
    }

    if (!igraph_is_directed(graph)) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (mode == IGRAPH_TO_UNDIRECTED_EACH) {
        igraph_es_t es;
        igraph_eit_t eit;

        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));
        IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
        IGRAPH_FINALLY(igraph_es_destroy, &es);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);

        while (!IGRAPH_EIT_END(eit)) {
            igraph_integer_t edge = IGRAPH_EIT_GET(eit);
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, IGRAPH_FROM(graph, edge)));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, IGRAPH_TO(graph, edge)));
            IGRAPH_EIT_NEXT(eit);
        }

        igraph_eit_destroy(&eit);
        igraph_es_destroy(&es);
        IGRAPH_FINALLY_CLEAN(2);

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   no_of_nodes,
                                   IGRAPH_UNDIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        igraph_vector_int_destroy(&edges);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, true, true, true);
        IGRAPH_FINALLY_CLEAN(2);
        igraph_destroy(graph);
        *graph = newgraph;

    } else if (mode == IGRAPH_TO_UNDIRECTED_COLLAPSE) {
        igraph_vector_int_t inadj, outadj;
        igraph_vector_int_t mergeinto;
        igraph_integer_t actedge = 0;

        if (attr) {
            IGRAPH_VECTOR_INT_INIT_FINALLY(&mergeinto, no_of_edges);
        }

        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));
        IGRAPH_VECTOR_INT_INIT_FINALLY(&inadj, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&outadj, 0);

        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_integer_t n_out, n_in;
            igraph_integer_t p1 = -1, p2 = -1;
            igraph_integer_t e1 = 0, e2 = 0, n1 = 0, n2 = 0, last;
            IGRAPH_CHECK(igraph_incident(graph, &outadj, i, IGRAPH_OUT));
            IGRAPH_CHECK(igraph_incident(graph, &inadj, i, IGRAPH_IN));
            n_out = igraph_vector_int_size(&outadj);
            n_in = igraph_vector_int_size(&inadj);

#define STEPOUT() if ( (++p1) < n_out) {    \
        e1 = VECTOR(outadj)[p1]; \
        n1 = IGRAPH_TO(graph, e1);      \
    }
#define STEPIN()  if ( (++p2) < n_in) {         \
        e2 = VECTOR(inadj )[p2]; \
        n2 = IGRAPH_FROM(graph, e2);        \
    }
#define ADD_NEW_EDGE() { \
    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i)); \
    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, last)); \
}
#define MERGE_INTO_CURRENT_EDGE(which) { \
    if (attr) { \
        VECTOR(mergeinto)[which] = actedge; \
    } \
}

            STEPOUT();
            STEPIN();

            while (p1 < n_out && n1 <= i && p2 < n_in && n2 <= i) {
                last = (n1 <= n2) ? n1 : n2;
                ADD_NEW_EDGE();
                while (p1 < n_out && last == n1) {
                    MERGE_INTO_CURRENT_EDGE(e1);
                    STEPOUT();
                }
                while (p2 < n_in && last == n2) {
                    MERGE_INTO_CURRENT_EDGE(e2);
                    STEPIN();
                }
                actedge++;
            }

            while (p1 < n_out && n1 <= i) {
                last = n1;
                ADD_NEW_EDGE();
                while (p1 < n_out && last == n1) {
                    MERGE_INTO_CURRENT_EDGE(e1);
                    STEPOUT();
                }
                actedge++;
            }

            while (p2 < n_in && n2 <= i) {
                last = n2;
                ADD_NEW_EDGE();
                while (p2 < n_in && last == n2) {
                    MERGE_INTO_CURRENT_EDGE(e2);
                    STEPIN();
                }
                actedge++;
            }
        }

#undef MERGE_INTO_CURRENT_EDGE
#undef ADD_NEW_EDGE
#undef STEPOUT
#undef STEPIN

        igraph_vector_int_destroy(&outadj);
        igraph_vector_int_destroy(&inadj);
        IGRAPH_FINALLY_CLEAN(2);

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   no_of_nodes,
                                   IGRAPH_UNDIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        igraph_vector_int_destroy(&edges);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, true, true, /*edges*/ false); /* no edge attributes */

        if (attr) {
            igraph_fixed_vectorlist_t vl;
            IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto, actedge));
            IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

            IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &newgraph, &vl.vecs, edge_comb));

            igraph_fixed_vectorlist_destroy(&vl);
            IGRAPH_FINALLY_CLEAN(1);
        }

        IGRAPH_FINALLY_CLEAN(2);
        igraph_destroy(graph);
        *graph = newgraph;

        if (attr) {
            igraph_vector_int_destroy(&mergeinto);
            IGRAPH_FINALLY_CLEAN(1);
        }
    } else if (mode == IGRAPH_TO_UNDIRECTED_MUTUAL) {
        igraph_vector_int_t inadj, outadj;
        igraph_vector_int_t mergeinto;
        igraph_integer_t actedge = 0;

        if (attr) {
            IGRAPH_VECTOR_INT_INIT_FINALLY(&mergeinto, no_of_edges);
            igraph_vector_int_fill(&mergeinto, -1);
        }

        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));
        IGRAPH_VECTOR_INT_INIT_FINALLY(&inadj, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&outadj, 0);

        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_integer_t n_out, n_in;
            igraph_integer_t p1 = -1, p2 = -1;
            igraph_integer_t e1 = 0, e2 = 0, n1 = 0, n2 = 0;
            IGRAPH_CHECK(igraph_incident(graph, &outadj, i,
                                         IGRAPH_OUT));
            IGRAPH_CHECK(igraph_incident(graph, &inadj,  i,
                                         IGRAPH_IN));
            n_out = igraph_vector_int_size(&outadj);
            n_in = igraph_vector_int_size(&inadj);

#define STEPOUT() if ( (++p1) < n_out) {    \
        e1 = VECTOR(outadj)[p1]; \
        n1 = IGRAPH_TO(graph, e1);      \
    }
#define STEPIN()  if ( (++p2) < n_in) {         \
        e2 = VECTOR(inadj )[p2]; \
        n2 = IGRAPH_FROM(graph, e2);        \
    }

            STEPOUT();
            STEPIN();

            while (p1 < n_out && n1 <= i && p2 < n_in && n2 <= i) {
                if (n1 == n2) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, n1));
                    if (attr) {
                        VECTOR(mergeinto)[e1] = actedge;
                        VECTOR(mergeinto)[e2] = actedge;
                        actedge++;
                    }
                    STEPOUT();
                    STEPIN();
                } else if (n1 < n2) {
                    STEPOUT();
                } else { /* n2<n1 */
                    STEPIN();
                }
            }
        }

#undef STEPOUT
#undef STEPIN

        igraph_vector_int_destroy(&outadj);
        igraph_vector_int_destroy(&inadj);
        IGRAPH_FINALLY_CLEAN(2);

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   no_of_nodes,
                                   IGRAPH_UNDIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        igraph_vector_int_destroy(&edges);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, true, true, /*edges*/ false); /* no edge attributes */

        if (attr) {
            igraph_fixed_vectorlist_t vl;
            IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto, actedge));
            IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

            IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &newgraph, &vl.vecs, edge_comb));

            igraph_fixed_vectorlist_destroy(&vl);
            IGRAPH_FINALLY_CLEAN(1);
        }

        IGRAPH_FINALLY_CLEAN(2);
        igraph_destroy(graph);
        *graph = newgraph;

        if (attr) {
            igraph_vector_int_destroy(&mergeinto);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return IGRAPH_SUCCESS;
}

#define WEIGHT_OF(eid) (weights ? VECTOR(*weights)[eid] : 1)

/**
 * \function igraph_get_stochastic
 * \brief Stochastic adjacency matrix of a graph.
 *
 * Stochastic matrix of a graph. The stochastic matrix of a graph is
 * its adjacency matrix, normalized row-wise or column-wise, such that
 * the sum of each row (or column) is one.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized matrix, the result is stored here.
 *   It will be resized as needed.
 * \param column_wise Whether to normalize column-wise.
 * \return Error code.
 *
 * Time complexity: O(|V||V|), |V| is the number of vertices in the graph.
 *
 * \sa \ref igraph_get_stochastic_sparse(), the sparse version of this
 * function.
 */

igraph_error_t igraph_get_stochastic(
    const igraph_t *graph, igraph_matrix_t *res, igraph_bool_t column_wise,
    const igraph_vector_t *weights
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t i, from, to;
    igraph_vector_t sums;
    igraph_real_t sum;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);

    IGRAPH_VECTOR_INIT_FINALLY(&sums, no_of_nodes);

    if (directed) {
        IGRAPH_CHECK(igraph_strength(
            graph, &sums, igraph_vss_all(),
            column_wise ? IGRAPH_IN : IGRAPH_OUT,
            /* loops = */ true, weights
        ));

        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            sum = VECTOR(sums)[column_wise ? to : from];
            MATRIX(*res, from, to) += WEIGHT_OF(i) / sum;
        }
    } else {
        IGRAPH_CHECK(igraph_strength(
            graph, &sums, igraph_vss_all(), IGRAPH_ALL,
            /* loops = */ true, weights
        ));

        for (i = 0; i < no_of_edges; i++) {
            from = IGRAPH_FROM(graph, i);
            to = IGRAPH_TO(graph, i);
            MATRIX(*res, from, to) += WEIGHT_OF(i) / VECTOR(sums)[column_wise ? to : from];
            MATRIX(*res, to, from) += WEIGHT_OF(i) / VECTOR(sums)[column_wise ? from: to];
        }
    }

    igraph_vector_destroy(&sums);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

#undef WEIGHT_OF

/**
 * \function igraph_get_stochastic_sparse
 * \brief The stochastic adjacency matrix of a graph.
 *
 * Stochastic matrix of a graph. The stochastic matrix of a graph is
 * its adjacency matrix, normalized row-wise or column-wise, such that
 * the sum of each row (or column) is one.
 *
 * \param graph The input graph.
 * \param res Pointer to an \em initialized sparse matrix, the
 *    result is stored here. The matrix will be resized as needed.
 * \param column_wise Whether to normalize column-wise.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa \ref igraph_get_stochastic(), the dense version of this function.
 */

igraph_error_t igraph_get_stochastic_sparse(
    const igraph_t *graph, igraph_sparsemat_t *res, igraph_bool_t column_wise,
    const igraph_vector_t *weights
) {
    IGRAPH_CHECK(igraph_get_adjacency_sparse(graph, res, IGRAPH_GET_ADJACENCY_BOTH, weights, IGRAPH_LOOPS_TWICE));

    if (column_wise) {
        IGRAPH_CHECK(igraph_sparsemat_normalize_cols(res, /* allow_zeros = */ false));
    } else {
        IGRAPH_CHECK(igraph_sparsemat_normalize_rows(res, /* allow_zeros = */ false));
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_get_stochastic_sparsemat
 * \brief Stochastic adjacency matrix of a graph (deprecated).
 *
 * This function is deprecated in favour of \ref igraph_get_stochastic_sparse(),
 * but does not work in an identical way. This function takes an \em uninitialized
 * \c igraph_sparsemat_t while \ref igraph_get_stochastic_sparse() takes
 * an already initialized one.
 *
 * \param graph The input graph.
 * \param res Pointer to an \em uninitialized sparse matrix, the
 *    result is stored here. The matrix will be resized as needed.
 * \param column_wise Whether to normalize column-wise. For undirected
 *    graphs this argument does not have any effect.
 * \return Error code.
 *
 * \deprecated-by igraph_get_stochastic_sparse 0.10.0
 */

igraph_error_t igraph_get_stochastic_sparsemat(const igraph_t *graph,
                                               igraph_sparsemat_t *res,
                                               igraph_bool_t column_wise) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t nzmax = igraph_is_directed(graph) ? no_of_edges : 2*no_of_edges;
    IGRAPH_CHECK(igraph_sparsemat_init(res, no_of_nodes, no_of_nodes, nzmax));
    return igraph_get_stochastic_sparse(graph, res, column_wise, NULL);
}


/**
 * \ingroup conversion
 * \function igraph_to_prufer
 * \brief Converts a tree to its Pr&uuml;fer sequence.
 *
 * A Pr&uuml;fer sequence is a unique sequence of integers associated
 * with a labelled tree. A tree on n >= 2 vertices can be represented by a
 * sequence of n-2 integers, each between 0 and n-1 (inclusive).
 *
 * \param graph Pointer to an initialized graph object which
          must be a tree on n >= 2 vertices.
 * \param prufer A pointer to the integer vector that should hold the Pr&uuml;fer sequence;
                 the vector must be initialized and will be resized to n - 2.
 * \return Error code:
 *          \clist
 *          \cli IGRAPH_ENOMEM
 *             there is not enough memory to perform the operation.
 *          \cli IGRAPH_EINVAL
 *             the graph is not a tree or it is has less than vertices
 *          \endclist
 *
 * \sa \ref igraph_from_prufer()
 *
 */
igraph_error_t igraph_to_prufer(const igraph_t *graph, igraph_vector_int_t* prufer) {
    /* For generating the Prüfer sequence, we enumerate the vertices u of the tree.
       We keep track of the degrees of all vertices, treating vertices
       of degree 0 as removed. We maintain the invariant that all leafs
       that are still contained in the tree are >= u.
       If u is a leaf, we remove it and add its unique neighbor to the Prüfer
       sequence. If the removal of u turns the neighbor into a leaf which is < u,
       we repeat the procedure for the new leaf and so on. */
    igraph_integer_t u;
    igraph_vector_int_t degrees;
    igraph_vector_int_t neighbors;
    igraph_integer_t prufer_index = 0;
    igraph_integer_t n = igraph_vcount(graph);
    igraph_bool_t is_tree = false;

    IGRAPH_CHECK(igraph_is_tree(graph, &is_tree, NULL, IGRAPH_ALL));

    if (!is_tree) {
        IGRAPH_ERROR("The graph must be a tree", IGRAPH_EINVAL);
    }

    if (n < 2) {
        IGRAPH_ERROR("The tree must have at least 2 vertices", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(prufer, n - 2));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, n);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neighbors, 1);

    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS));

    for (u = 0; u < n; ++u) {
        igraph_integer_t degree = VECTOR(degrees)[u];
        igraph_integer_t leaf = u;

        while (degree == 1 && leaf <= u) {
            igraph_integer_t neighbor = 0;
            igraph_integer_t neighbor_count = 0;

            VECTOR(degrees)[leaf] = 0; /* mark leaf v as deleted */

            IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, leaf, IGRAPH_ALL));

            /* Find the unique remaining neighbor of the leaf */
            neighbor_count = igraph_vector_int_size(&neighbors);
            for (igraph_integer_t i = 0; i < neighbor_count; i++) {
                neighbor = VECTOR(neighbors)[i];
                if (VECTOR(degrees)[neighbor] > 0) {
                    break;
                }
            }

            /* remember that we have removed the leaf */
            VECTOR(degrees)[neighbor]--;
            degree = VECTOR(degrees)[neighbor];

            /* Add the neighbor to the prufer sequence unless it is the last vertex
            (i.e. degree == 0) */
            if (degree > 0) {
                VECTOR(*prufer)[prufer_index] = neighbor;
                prufer_index++;
            }
            leaf = neighbor;
        }
    }

    igraph_vector_int_destroy(&degrees);
    igraph_vector_int_destroy(&neighbors);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
