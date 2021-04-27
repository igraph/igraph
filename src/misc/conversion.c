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
#include "misc/conversion_internal.h"

/**
 * \ingroup conversion
 * \function igraph_get_adjacency
 * \brief Returns the adjacency matrix of a graph
 *
 * </para><para>
 * The result is an adjacency matrix. Entry i, j of the matrix
 * contains the number of edges connecting vertex i to vertex j.
 * \param graph Pointer to the graph to convert
 * \param res Pointer to an initialized matrix object, it will be
 *        resized if needed.
 * \param type Constant giving the type of the adjacency matrix to
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
 * \param type eids Logical, if true, then the edges ids plus one
 *        are stored in the adjacency matrix, instead of the number of
 *        edges between the two vertices. (The plus one is needed, since
 *        edge ids start from zero, and zero means no edge in this case.)
 * \return Error code:
 *        \c IGRAPH_EINVAL invalid type argument.
 *
 * \sa igraph_get_adjacency_sparse if you want a sparse matrix representation
 *
 * Time complexity: O(|V||V|),
 * |V| is the
 * number of vertices in the graph.
 */

int igraph_get_adjacency(const igraph_t *graph, igraph_matrix_t *res,
                         igraph_get_adjacency_t type, igraph_bool_t eids) {

    igraph_eit_t edgeit;
    long int no_of_nodes = igraph_vcount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    int retval = 0;
    long int from, to;
    igraph_integer_t ffrom, fto;

    IGRAPH_CHECK(igraph_matrix_resize(res, no_of_nodes, no_of_nodes));
    igraph_matrix_null(res);
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
    IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);

    if (directed) {
        while (!IGRAPH_EIT_END(edgeit)) {
            long int edge = IGRAPH_EIT_GET(edgeit);
            igraph_edge(graph, (igraph_integer_t) edge, &ffrom, &fto);
            from = ffrom;
            to = fto;
            if (eids) {
                MATRIX(*res, from, to) = edge + 1;
            } else {
                MATRIX(*res, from, to) += 1;
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else if (type == IGRAPH_GET_ADJACENCY_UPPER) {
        while (!IGRAPH_EIT_END(edgeit)) {
            long int edge = IGRAPH_EIT_GET(edgeit);
            igraph_edge(graph, (igraph_integer_t) edge, &ffrom, &fto);
            from = ffrom;
            to = fto;
            if (to < from) {
                if (eids) {
                    MATRIX(*res, to, from) = edge + 1;
                } else {
                    MATRIX(*res, to, from) += 1;
                }
            } else {
                if (eids) {
                    MATRIX(*res, from, to) = edge + 1;
                } else {
                    MATRIX(*res, from, to) += 1;
                }
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else if (type == IGRAPH_GET_ADJACENCY_LOWER) {
        while (!IGRAPH_EIT_END(edgeit)) {
            long int edge = IGRAPH_EIT_GET(edgeit);
            igraph_edge(graph, (igraph_integer_t) edge, &ffrom, &fto);
            from = ffrom;
            to = fto;
            if (to < from) {
                if (eids) {
                    MATRIX(*res, from, to) = edge + 1;
                } else {
                    MATRIX(*res, from, to) += 1;
                }
            } else {
                if (eids) {
                    MATRIX(*res, to, from) = edge + 1;
                } else {
                    MATRIX(*res, to, from) += 1;
                }
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else if (type == IGRAPH_GET_ADJACENCY_BOTH) {
        while (!IGRAPH_EIT_END(edgeit)) {
            long int edge = IGRAPH_EIT_GET(edgeit);
            igraph_edge(graph, (igraph_integer_t) edge, &ffrom, &fto);
            from = ffrom;
            to = fto;
            if (eids) {
                MATRIX(*res, from, to) = edge + 1;
            } else {
                MATRIX(*res, from, to) += 1;
            }
            if (from != to) {
                if (eids) {
                    MATRIX(*res, to, from) = edge + 1;
                } else {
                    MATRIX(*res, to, from) += 1;
                }
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else {
        IGRAPH_ERROR("Invalid type argument", IGRAPH_EINVAL);
    }

    igraph_eit_destroy(&edgeit);
    IGRAPH_FINALLY_CLEAN(1);
    return retval;
}

/**
 * \ingroup conversion
 * \function igraph_get_adjacency_sparse
 * \brief Returns the adjacency matrix of a graph in sparse matrix format.
 *
 * </para><para>
 * The result is an adjacency matrix. Entry i, j of the matrix
 * contains the number of edges connecting vertex i to vertex j.
 * \param graph Pointer to the graph to convert.
 * \param res Pointer to an initialized sparse matrix object, it will be
 *        resized if needed.
 * \param type Constant giving the type of the adjacency matrix to
 *        create for undirected graphs. It is ignored for directed
 *        graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_GET_ADJACENCY_UPPER
 *          the upper right triangle of the matrix is used.
 *        \cli IGRAPH_GET_ADJACENCY_LOWER
 *          the lower left triangle of the matrix is used.
 *        \cli IGRAPH_GET_ADJACENCY_BOTH
 *          the whole matrix is used, a symmetric matrix is returned.
 *        \endclist
 * \return Error code:
 *        \c IGRAPH_EINVAL invalid type argument.
 *
 * \sa igraph_get_adjacency if you would like to get a normal matrix
 *   ( \type igraph_matrix_t )
 *
 * Time complexity: O(|V||V|),
 * |V| is the
 * number of vertices in the graph.
 */

int igraph_get_adjacency_sparse(const igraph_t *graph, igraph_spmatrix_t *res,
                                igraph_get_adjacency_t type) {

    igraph_eit_t edgeit;
    long int no_of_nodes = igraph_vcount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    long int from, to;
    igraph_integer_t ffrom, fto;

    igraph_spmatrix_null(res);
    IGRAPH_CHECK(igraph_spmatrix_resize(res, no_of_nodes, no_of_nodes));
    IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &edgeit));
    IGRAPH_FINALLY(igraph_eit_destroy, &edgeit);

    if (directed) {
        while (!IGRAPH_EIT_END(edgeit)) {
            igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
            from = ffrom;
            to = fto;
            igraph_spmatrix_add_e(res, from, to, 1);
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else if (type == IGRAPH_GET_ADJACENCY_UPPER) {
        while (!IGRAPH_EIT_END(edgeit)) {
            igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
            from = ffrom;
            to = fto;
            if (to < from) {
                igraph_spmatrix_add_e(res, to, from, 1);
            } else {
                igraph_spmatrix_add_e(res, from, to, 1);
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else if (type == IGRAPH_GET_ADJACENCY_LOWER) {
        while (!IGRAPH_EIT_END(edgeit)) {
            igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
            from = ffrom;
            to = fto;
            if (to > from) {
                igraph_spmatrix_add_e(res, to, from, 1);
            } else {
                igraph_spmatrix_add_e(res, from, to, 1);
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else if (type == IGRAPH_GET_ADJACENCY_BOTH) {
        while (!IGRAPH_EIT_END(edgeit)) {
            igraph_edge(graph, IGRAPH_EIT_GET(edgeit), &ffrom, &fto);
            from = ffrom;
            to = fto;
            igraph_spmatrix_add_e(res, from, to, 1);
            if (from != to) {
                igraph_spmatrix_add_e(res, to, from, 1);
            }
            IGRAPH_EIT_NEXT(edgeit);
        }
    } else {
        IGRAPH_ERROR("Invalid type argument.", IGRAPH_EINVAL);
    }

    igraph_eit_destroy(&edgeit);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup conversion
 * \function igraph_get_edgelist
 * \brief Returns the list of edges in a graph
 *
 * </para><para>The order of the edges is given by the edge ids.
 * \param graph Pointer to the graph object
 * \param res Pointer to an initialized vector object, it will be
 *        resized.
 * \param bycol Logical, if true, the edges will be returned
 *        columnwise, e.g. the first edge is
 *        <code>res[0]->res[|E|]</code>, the second is
 *        <code>res[1]->res[|E|+1]</code>, etc.
 * \return Error code.
 *
 * \sa \ref igraph_edges() to return the result only for some edge ids.
 *
 * Time complexity: O(|E|), the
 * number of edges in the graph.
 */

int igraph_get_edgelist(const igraph_t *graph, igraph_vector_t *res, igraph_bool_t bycol) {

    igraph_eit_t edgeit;
    long int no_of_edges = igraph_ecount(graph);
    long int vptr = 0;
    igraph_integer_t from, to;

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_edges * 2));
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
    return 0;
}

/**
 * \function igraph_to_directed
 * \brief Convert an undirected graph to a directed one
 *
 * </para><para>
 * If the supplied graph is directed, this function does nothing.
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

int igraph_to_directed(igraph_t *graph,
                       igraph_to_directed_t mode) {
    long int no_of_edges = igraph_ecount(graph);
    long int no_of_nodes = igraph_vcount(graph);

    if (igraph_is_directed(graph)) {
        return IGRAPH_SUCCESS;
    }   

    switch (mode) {
    case IGRAPH_TO_DIRECTED_ARBITRARY:
    case IGRAPH_TO_DIRECTED_RANDOM:
    case IGRAPH_TO_DIRECTED_ACYCLIC:
      {
        igraph_t newgraph;
        igraph_vector_t edges;
        long int size = no_of_edges * 2;
        long int i;

        IGRAPH_VECTOR_INIT_FINALLY(&edges, size);
        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));

        if (mode == IGRAPH_TO_DIRECTED_RANDOM) {
            RNG_BEGIN();

            for (i=0; i < no_of_edges; ++i) {
                if (RNG_INTEGER(0,1)) {
                    igraph_real_t temp = VECTOR(edges)[2*i];
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
            for (i=0; i < no_of_edges; ++i) {
                if (VECTOR(edges)[2*i] > VECTOR(edges)[2*i+1]) {
                    igraph_real_t temp = VECTOR(edges)[2*i];
                    VECTOR(edges)[2*i] = VECTOR(edges)[2*i+1];
                    VECTOR(edges)[2*i+1] = temp;
                }
            }
        }

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   (igraph_integer_t) no_of_nodes,
                                   IGRAPH_DIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);        
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 1);
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(2);

        igraph_destroy(graph);
        *graph = newgraph;

        break;
      }
    case IGRAPH_TO_DIRECTED_MUTUAL:
      {
        igraph_t newgraph;
        igraph_vector_t edges;
        igraph_vector_t index;
        long int size = no_of_edges * 4;
        long int i;
        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&edges, size));
        IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, 0));
        IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges * 4));
        IGRAPH_VECTOR_INIT_FINALLY(&index, no_of_edges * 2);
        for (i = 0; i < no_of_edges; i++) {
            VECTOR(edges)[no_of_edges * 2 + i * 2]  = VECTOR(edges)[i * 2 + 1];
            VECTOR(edges)[no_of_edges * 2 + i * 2 + 1] = VECTOR(edges)[i * 2];
            VECTOR(index)[i] = VECTOR(index)[no_of_edges + i] = i;
        }

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   (igraph_integer_t) no_of_nodes,
                                   IGRAPH_DIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1,/*edges=*/0);
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, &newgraph, &index));

        igraph_vector_destroy(&index);
        igraph_vector_destroy(&edges);        
        IGRAPH_FINALLY_CLEAN(3);

        igraph_destroy(graph);
        *graph = newgraph;

        break;
      }
    default:
        IGRAPH_ERROR("Cannot direct graph, invalid mode", IGRAPH_EINVAL);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_to_undirected
 * \brief Convert a directed graph to an undirected one.
 *
 * </para><para>
 * If the supplied graph is undirected, this function does nothing.
 * \param graph The graph object to convert.
 * \param mode Constant, specifies the details of how exactly the
 *        conversion is done. Possible values: \c
 *        IGRAPH_TO_UNDIRECTED_EACH: the number of edges remains
 *        constant, an undirected edge is created for each directed
 *        one, this version might create graphs with multiple edges;
 *        \c IGRAPH_TO_UNDIRECTED_COLLAPSE: one undirected edge will
 *        be created for each pair of vertices which are connected
 *        with at least one directed edge, no multiple edges will be
 *        created. \c IGRAPH_TO_UNDIRECTED_MUTUAL creates an undirected
 *        edge for each pair of mutual edges in the directed graph.
 *        Non-mutual edges are lost. This mode might create multiple
 *        edges.
 * \param edge_comb What to do with the edge attributes. See the igraph
 *        manual section about attributes for details.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 *
 * \example examples/simple/igraph_to_undirected.c
 */

int igraph_to_undirected(igraph_t *graph,
                         igraph_to_undirected_t mode,
                         const igraph_attribute_combination_t *edge_comb) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vector_t edges;
    igraph_t newgraph;
    igraph_bool_t attr = edge_comb && igraph_has_attribute_table();

    if (mode != IGRAPH_TO_UNDIRECTED_EACH &&
        mode != IGRAPH_TO_UNDIRECTED_COLLAPSE &&
        mode != IGRAPH_TO_UNDIRECTED_MUTUAL) {
        IGRAPH_ERROR("Cannot undirect graph, invalid mode", IGRAPH_EINVAL);
    }

    if (!igraph_is_directed(graph)) {
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    if (mode == IGRAPH_TO_UNDIRECTED_EACH) {
        igraph_es_t es;
        igraph_eit_t eit;

        IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));
        IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
        IGRAPH_FINALLY(igraph_es_destroy, &es);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);

        while (!IGRAPH_EIT_END(eit)) {
            long int edge = IGRAPH_EIT_GET(eit);
            igraph_integer_t from, to;
            igraph_edge(graph, (igraph_integer_t) edge, &from, &to);
            IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
            IGRAPH_EIT_NEXT(eit);
        }

        igraph_eit_destroy(&eit);
        igraph_es_destroy(&es);
        IGRAPH_FINALLY_CLEAN(2);

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   (igraph_integer_t) no_of_nodes,
                                   IGRAPH_UNDIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        igraph_vector_destroy(&edges);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 1);
        IGRAPH_FINALLY_CLEAN(2);
        igraph_destroy(graph);
        *graph = newgraph;

    } else if (mode == IGRAPH_TO_UNDIRECTED_COLLAPSE) {
        igraph_vector_t inadj, outadj;
        long int i;
        igraph_vector_t mergeinto;
        long int actedge = 0;

        if (attr) {
            IGRAPH_VECTOR_INIT_FINALLY(&mergeinto, no_of_edges);
        }

        IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));
        IGRAPH_VECTOR_INIT_FINALLY(&inadj, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&outadj, 0);

        for (i = 0; i < no_of_nodes; i++) {
            long int n_out, n_in;
            long int p1 = -1, p2 = -1;
            long int e1 = 0, e2 = 0, n1 = 0, n2 = 0;
            IGRAPH_CHECK(igraph_incident(graph, &outadj, (igraph_integer_t) i,
                                         IGRAPH_OUT));
            IGRAPH_CHECK(igraph_incident(graph, &inadj, (igraph_integer_t) i,
                                         IGRAPH_IN));
            n_out = igraph_vector_size(&outadj);
            n_in = igraph_vector_size(&inadj);

#define STEPOUT() if ( (++p1) < n_out) {    \
        e1 = (long int) VECTOR(outadj)[p1]; \
        n1 = IGRAPH_TO(graph, e1);      \
    }
#define STEPIN()  if ( (++p2) < n_in) {         \
        e2 = (long int) VECTOR(inadj )[p2]; \
        n2 = IGRAPH_FROM(graph, e2);        \
    }

            STEPOUT();
            STEPIN();

            while (p1 < n_out && n1 <= i && p2 < n_in && n2 <= i) {
                long int last;
                if (n1 == n2) {
                    last = n1;
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, n1));
                    if (attr) {
                        VECTOR(mergeinto)[e1] = actedge;
                        VECTOR(mergeinto)[e2] = actedge;
                        actedge++;
                    }
                    while (p1 < n_out && last == n1) {
                        STEPOUT();
                    }
                    while (p2 < n_in  && last == n2) {
                        STEPIN ();
                    }
                } else if (n1 < n2) {
                    last = n1;
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, n1));
                    if (attr) {
                        VECTOR(mergeinto)[e1] = actedge;
                        actedge++;
                    }
                    while (p1 < n_out && last == n1) {
                        STEPOUT();
                    }
                } else { /* n2<n1 */
                    last = n2;
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, n2));
                    if (attr) {
                        VECTOR(mergeinto)[e2] = actedge;
                        actedge++;
                    }
                    while (p2 < n_in && last == n2) {
                        STEPIN();
                    }
                }
            }
            while (p1 < n_out && n1 <= i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, n1));
                if (attr) {
                    VECTOR(mergeinto)[e1] = actedge;
                    actedge++;
                }
                STEPOUT();
            }
            while (p2 < n_in && n2 <= i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, n2));
                if (attr) {
                    VECTOR(mergeinto)[e2] = actedge;
                    actedge++;
                }
                STEPIN();
            }
        }

#undef STEPOUT
#undef STEPIN

        igraph_vector_destroy(&outadj);
        igraph_vector_destroy(&inadj);
        IGRAPH_FINALLY_CLEAN(2);

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   (igraph_integer_t) no_of_nodes,
                                   IGRAPH_UNDIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        igraph_vector_destroy(&edges);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 0); /* no edge attributes */

        if (attr) {
            igraph_fixed_vectorlist_t vl;
            IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto,
                         actedge));
            IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

            IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &newgraph, &vl.v,
                         edge_comb));

            igraph_fixed_vectorlist_destroy(&vl);
            IGRAPH_FINALLY_CLEAN(1);
        }

        IGRAPH_FINALLY_CLEAN(2);
        igraph_destroy(graph);
        *graph = newgraph;

        if (attr) {
            igraph_vector_destroy(&mergeinto);
            IGRAPH_FINALLY_CLEAN(1);
        }
    } else if (mode == IGRAPH_TO_UNDIRECTED_MUTUAL) {
        igraph_vector_t inadj, outadj;
        long int i;
        igraph_vector_t mergeinto;
        long int actedge = 0;

        if (attr) {
            IGRAPH_VECTOR_INIT_FINALLY(&mergeinto, no_of_edges);
            igraph_vector_fill(&mergeinto, -1);
        }

        IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));
        IGRAPH_VECTOR_INIT_FINALLY(&inadj, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&outadj, 0);

        for (i = 0; i < no_of_nodes; i++) {
            long int n_out, n_in;
            long int p1 = -1, p2 = -1;
            long int e1 = 0, e2 = 0, n1 = 0, n2 = 0;
            IGRAPH_CHECK(igraph_incident(graph, &outadj, (igraph_integer_t) i,
                                         IGRAPH_OUT));
            IGRAPH_CHECK(igraph_incident(graph, &inadj,  (igraph_integer_t) i,
                                         IGRAPH_IN));
            n_out = igraph_vector_size(&outadj);
            n_in = igraph_vector_size(&inadj);

#define STEPOUT() if ( (++p1) < n_out) {    \
        e1 = (long int) VECTOR(outadj)[p1]; \
        n1 = IGRAPH_TO(graph, e1);      \
    }
#define STEPIN()  if ( (++p2) < n_in) {         \
        e2 = (long int) VECTOR(inadj )[p2]; \
        n2 = IGRAPH_FROM(graph, e2);        \
    }

            STEPOUT();
            STEPIN();

            while (p1 < n_out && n1 <= i && p2 < n_in && n2 <= i) {
                if (n1 == n2) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, n1));
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

        igraph_vector_destroy(&outadj);
        igraph_vector_destroy(&inadj);
        IGRAPH_FINALLY_CLEAN(2);

        IGRAPH_CHECK(igraph_create(&newgraph, &edges,
                                   (igraph_integer_t) no_of_nodes,
                                   IGRAPH_UNDIRECTED));
        IGRAPH_FINALLY(igraph_destroy, &newgraph);
        igraph_vector_destroy(&edges);
        IGRAPH_I_ATTRIBUTE_DESTROY(&newgraph);
        IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 1, 1, 0); /* no edge attributes */

        if (attr) {
            igraph_fixed_vectorlist_t vl;
            IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto,
                         actedge));
            IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

            IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &newgraph, &vl.v,
                         edge_comb));

            igraph_fixed_vectorlist_destroy(&vl);
            IGRAPH_FINALLY_CLEAN(1);
        }

        IGRAPH_FINALLY_CLEAN(2);
        igraph_destroy(graph);
        *graph = newgraph;

        if (attr) {
            igraph_vector_destroy(&mergeinto);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return 0;
}

/**
 * \function igraph_get_stochastic
 * Stochastic adjacency matrix of a graph
 *
 * Stochastic matrix of a graph. The stochastic matrix of a graph is
 * its adjacency matrix, normalized row-wise or column-wise, such that
 * the sum of each row (or column) is one.
 * \param graph The input graph.
 * \param sparsemat Pointer to an initialized matrix, the
 *    result is stored here.
 * \param column_wise Whether to normalize column-wise. For undirected
 *    graphs this argument does not have any effect.
 * \return Error code.
 *
 * Time complexity: O(|V||V|), quadratic in the number of vertices.
 *
 * \sa igraph_get_stochastic_sparsemat(), the sparse version of this
 * function.
 */

int igraph_get_stochastic(const igraph_t *graph,
                          igraph_matrix_t *matrix,
                          igraph_bool_t column_wise) {

    int no_of_nodes = igraph_vcount(graph);
    igraph_real_t sum;
    int i, j;

    IGRAPH_CHECK(igraph_get_adjacency(graph, matrix,
                                      IGRAPH_GET_ADJACENCY_BOTH, /*eids=*/ 0));

    if (!column_wise) {
        for (i = 0; i < no_of_nodes; i++) {
            sum = 0.0;
            for (j = 0; j < no_of_nodes; j++) {
                sum += MATRIX(*matrix, i, j);
            }
            for (j = 0; j < no_of_nodes; j++) {
                MATRIX(*matrix, i, j) /= sum;
            }
        }
    } else {
        for (i = 0; i < no_of_nodes; i++) {
            sum = 0.0;
            for (j = 0; j < no_of_nodes; j++) {
                sum += MATRIX(*matrix, j, i);
            }
            for (j = 0; j < no_of_nodes; j++) {
                MATRIX(*matrix, j, i) /= sum;
            }
        }
    }

    return 0;
}


int igraph_i_normalize_sparsemat(igraph_sparsemat_t *sparsemat,
                                 igraph_bool_t column_wise) {
    igraph_vector_t sum;
    int no_of_nodes = (int) igraph_sparsemat_nrow(sparsemat);
    int i;

    IGRAPH_VECTOR_INIT_FINALLY(&sum, no_of_nodes);

    if (!column_wise) {
        IGRAPH_CHECK(igraph_sparsemat_rowsums(sparsemat, &sum));
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(sum)[i] == 0.0) {
                IGRAPH_ERROR("Zero out-degree vertices not allowed",
                             IGRAPH_EINVAL);
            }
            VECTOR(sum)[i] = 1.0 / VECTOR(sum)[i];
        }
        IGRAPH_CHECK(igraph_sparsemat_scale_rows(sparsemat, &sum));
    } else {
        IGRAPH_CHECK(igraph_sparsemat_colsums(sparsemat, &sum));
        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(sum)[i] == 0.0) {
                IGRAPH_ERROR("Zero out-degree vertices not allowed",
                             IGRAPH_EINVAL);
            }
            VECTOR(sum)[i] = 1.0 / VECTOR(sum)[i];
        }
        IGRAPH_CHECK(igraph_sparsemat_scale_cols(sparsemat, &sum));
    }

    igraph_vector_destroy(&sum);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_get_stochastic_sparsemat
 * \brief Stochastic adjacency matrix of a graph
 *
 * Stochastic matrix of a graph. The stochastic matrix of a graph is
 * its adjacency matrix, normalized row-wise or column-wise, such that
 * the sum of each row (or column) is one.
 * \param graph The input graph.
 * \param sparsemat Pointer to an uninitialized sparse matrix, the
 *    result is stored here.
 * \param column_wise Whether to normalize column-wise. For undirected
 *    graphs this argument does not have any effect.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa igraph_get_stochastic(), the dense version of this function.
 */

int igraph_get_stochastic_sparsemat(const igraph_t *graph,
                                    igraph_sparsemat_t *sparsemat,
                                    igraph_bool_t column_wise) {

    IGRAPH_CHECK(igraph_get_sparsemat(graph, sparsemat));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, sparsemat);
    IGRAPH_CHECK(igraph_i_normalize_sparsemat(sparsemat, column_wise));
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


/**
 * \ingroup conversion
 * \function igraph_to_prufer
 * \brief Converts a tree to its Pr&uuml;fer sequence
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
int igraph_to_prufer(const igraph_t *graph, igraph_vector_int_t* prufer) {
    /* For generating the Prüfer sequence, we enumerate the vertices u of the tree.
       We keep track of the degrees of all vertices, treating vertices
       of degree 0 as removed. We maintain the invariant that all leafs
       that are still contained in the tree are >= u.
       If u is a leaf, we remove it and add its unique neighbor to the prüfer
       sequence. If the removal of u turns the neighbor into a leaf which is < u,
       we repeat the procedure for the new leaf and so on. */
    igraph_integer_t u;
    igraph_vector_t degrees, neighbors;
    igraph_integer_t prufer_index = 0;
    igraph_integer_t n = igraph_vcount(graph);
    igraph_bool_t is_tree = 0;

    IGRAPH_CHECK(igraph_is_tree(graph, &is_tree, NULL, IGRAPH_ALL));

    if (!is_tree) {
        IGRAPH_ERROR("The graph must be a tree", IGRAPH_EINVAL);
    }

    if (n < 2) {
        IGRAPH_ERROR("The tree must have at least 2 vertices", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_vector_int_resize(prufer, n - 2));
    IGRAPH_VECTOR_INIT_FINALLY(&degrees, n);
    IGRAPH_VECTOR_INIT_FINALLY(&neighbors, 1);

    IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS));

    for (u = 0; u < n; ++u) {
        igraph_integer_t degree = VECTOR(degrees)[u];
        igraph_integer_t leaf = u;

        while (degree == 1 && leaf <= u) {
            igraph_integer_t i;
            igraph_integer_t neighbor = 0;
            igraph_integer_t neighbor_count = 0;

            VECTOR(degrees)[leaf] = 0; /* mark leaf v as deleted */

            IGRAPH_CHECK(igraph_neighbors(graph, &neighbors, leaf, IGRAPH_ALL));

            /* Find the unique remaining neighbor of the leaf */
            neighbor_count = igraph_vector_size(&neighbors);
            for (i = 0; i < neighbor_count; i++) {
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

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&neighbors);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
