/*
   igraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_bipartite.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_constructors.h"
#include "igraph_dqueue.h"
#include "igraph_random.h"

#include "core/interruption.h"
#include "graph/attributes.h"
#include "internal/utils.h"
#include "math/safe_intop.h"
#include "misc/graphicality.h"
#include "random/random_internal.h"

/**
 * \section about_bipartite Bipartite networks in igraph
 *
 * <para>
 * A bipartite network contains two kinds of vertices and connections
 * are only possible between two vertices of different kinds. There are
 * many natural examples, e.g. movies and actors as vertices and a
 * movie is connected to all participating actors, etc.
 *
 * </para><para>
 * igraph does not have direct support for bipartite networks, at
 * least not at the C language level. In other words the igraph_t
 * structure does not contain information about the vertex types.
 * The C functions for bipartite networks usually have an additional
 * input argument to graph, called \c types, a boolean vector giving
 * the vertex types.
 *
 * </para><para>
 * Most functions creating bipartite networks are able to create this
 * extra vector, you just need to supply an initialized boolean vector
 * to them.</para>
 */

/**
 * \function igraph_bipartite_projection_size
 * \brief Calculate the number of vertices and edges in the bipartite projections.
 *
 * This function calculates the number of vertices and edges in the
 * two projections of a bipartite network. This is useful if you have
 * a big bipartite network and you want to estimate the amount of
 * memory you would need to calculate the projections themselves.
 *
 * \param graph The input graph.
 * \param types Boolean vector giving the vertex types of the graph.
 * \param vcount1 Pointer to an \c igraph_int_t, the number of
 *     vertices in the first projection is stored here. May be \c NULL
 *     if not needed.
 * \param ecount1 Pointer to an \c igraph_int_t, the number of
 *     edges in the first projection is stored here. May be \c NULL
 *     if not needed.
 * \param vcount2 Pointer to an \c igraph_int_t, the number of
 *     vertices in the second projection is stored here. May be \c NULL
 *     if not needed.
 * \param ecount2 Pointer to an \c igraph_int_t, the number of
 *     edges in the second projection is stored here. May be \c NULL
 *     if not needed.
 * \return Error code.
 *
 * \sa \ref igraph_bipartite_projection() to calculate the actual
 * projection.
 *
 * Time complexity: O(|V|*d^2+|E|), |V| is the number of vertices, |E|
 * is the number of edges, d is the average (total) degree of the
 * graphs.
 */

igraph_error_t igraph_bipartite_projection_size(const igraph_t *graph,
                                     const igraph_vector_bool_t *types,
                                     igraph_int_t *vcount1,
                                     igraph_int_t *ecount1,
                                     igraph_int_t *vcount2,
                                     igraph_int_t *ecount2) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t vc1 = 0, ec1 = 0, vc2 = 0, ec2 = 0;
    igraph_adjlist_t adjlist;
    igraph_vector_int_t added;

    if (igraph_vector_bool_size(types) != no_of_nodes) {
        IGRAPH_ERROR("Invalid bipartite type vector length.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&added, no_of_nodes);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *neis1;
        igraph_int_t neilen1;
        igraph_int_t *ecptr;
        if (VECTOR(*types)[i]) {
            vc2++;
            ecptr = &ec2;
        } else {
            vc1++;
            ecptr = &ec1;
        }
        neis1 = igraph_adjlist_get(&adjlist, i);
        neilen1 = igraph_vector_int_size(neis1);
        for (igraph_int_t j = 0; j < neilen1; j++) {
            igraph_int_t neilen2, nei = VECTOR(*neis1)[j];
            const igraph_vector_int_t *neis2 = igraph_adjlist_get(&adjlist, nei);
            if (IGRAPH_UNLIKELY(VECTOR(*types)[i] == VECTOR(*types)[nei])) {
                IGRAPH_ERROR("Non-bipartite edge found in bipartite projection.",
                             IGRAPH_EINVAL);
            }
            neilen2 = igraph_vector_int_size(neis2);
            for (igraph_int_t k = 0; k < neilen2; k++) {
                igraph_int_t nei2 = VECTOR(*neis2)[k];
                if (nei2 <= i) {
                    continue;
                }
                if (VECTOR(added)[nei2] == i + 1) {
                    continue;
                }
                VECTOR(added)[nei2] = i + 1;
                (*ecptr)++;
            }
        }
    }

    if (vcount1) {
        *vcount1 = vc1;
    }

    if (ecount1) {
        *ecount1 = ec1;
    }

    if (vcount2) {
        *vcount2 = vc2;
    }

    if (ecount2) {
        *ecount2 = ec2;
    }

    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&added);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_bipartite_projection(const igraph_t *graph,
                                         const igraph_vector_bool_t *types,
                                         igraph_t *proj,
                                         int which,
                                         igraph_vector_int_t *multiplicity) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t remaining_nodes = 0;
    igraph_vector_int_t vertex_perm, vertex_index;
    igraph_vector_int_t edges;
    igraph_adjlist_t adjlist;
    const igraph_vector_int_t *neis1, *neis2;
    igraph_int_t neilen1, neilen2;
    igraph_vector_int_t added;
    igraph_vector_int_t mult;

    if (which < 0) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_perm, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&vertex_perm, no_of_nodes));

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vertex_index, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&added, no_of_nodes);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* we won't need the 'mult' vector if 'multiplicity' is NULL, but MSVC will
     * throw warnings in the compiler output if we initialize it conditionally */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&mult, multiplicity ? no_of_nodes : 1);
    if (multiplicity) {
        igraph_vector_int_clear(multiplicity);
    }

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] == which) {
            VECTOR(vertex_index)[i] = remaining_nodes++;
            igraph_vector_int_push_back(&vertex_perm, i);
        }
    }

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] == which) {
            igraph_int_t new_i = VECTOR(vertex_index)[i];
            igraph_int_t iedges = 0;
            neis1 = igraph_adjlist_get(&adjlist, i);
            neilen1 = igraph_vector_int_size(neis1);
            for (igraph_int_t j = 0; j < neilen1; j++) {
                igraph_int_t nei = VECTOR(*neis1)[j];
                if (IGRAPH_UNLIKELY(VECTOR(*types)[i] == VECTOR(*types)[nei])) {
                    IGRAPH_ERROR("Non-bipartite edge found in bipartite projection.", IGRAPH_EINVAL);
                }
                neis2 = igraph_adjlist_get(&adjlist, nei);
                neilen2 = igraph_vector_int_size(neis2);
                for (igraph_int_t k = 0; k < neilen2; k++) {
                    igraph_int_t nei2 = VECTOR(*neis2)[k], new_nei2;
                    if (nei2 <= i) {
                        continue;
                    }
                    if (VECTOR(added)[nei2] == i + 1) {
                        if (multiplicity) {
                            VECTOR(mult)[nei2] += 1;
                        }
                        continue;
                    }
                    VECTOR(added)[nei2] = i + 1;
                    if (multiplicity) {
                        VECTOR(mult)[nei2] = 1;
                    }
                    iedges++;

                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, new_i));
                    if (multiplicity) {
                        /* If we need the multiplicity as well, then we put in the
                           old vertex IDs here and rewrite it later */
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, nei2));
                    } else {
                        new_nei2 = VECTOR(vertex_index)[nei2];
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, new_nei2));
                    }
                }
            }
            if (multiplicity) {
                /* OK, we need to go through all the edges added for vertex new_i
                   and check their multiplicity */
                igraph_int_t now = igraph_vector_int_size(&edges);
                igraph_int_t from = now - iedges * 2;
                for (igraph_int_t j = from; j < now; j += 2) {
                    igraph_int_t nei2 = VECTOR(edges)[j + 1];
                    igraph_int_t new_nei2 = VECTOR(vertex_index)[nei2];
                    igraph_int_t m = VECTOR(mult)[nei2];
                    VECTOR(edges)[j + 1] = new_nei2;
                    IGRAPH_CHECK(igraph_vector_int_push_back(multiplicity, m));
                }
            }
        } /* if VECTOR(*type)[i] == which */
    }

    igraph_vector_int_destroy(&mult);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&added);
    igraph_vector_int_destroy(&vertex_index);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(proj, &edges, remaining_nodes, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, proj);

    /* copy graph attributes */
    IGRAPH_CHECK(igraph_i_attribute_copy(proj, graph, true, /* vertex= */ false, /* edge= */ false));

    /* copy vertex attributes */
    IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, proj, &vertex_perm));

    igraph_vector_int_destroy(&vertex_perm);
    IGRAPH_FINALLY_CLEAN(2); /* +1 for proj1 */

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bipartite_projection
 * \brief Create one or both projections of a bipartite (two-mode) network.
 *
 * Creates one or both projections of a bipartite graph.
 *
 * </para><para>
 * A graph is called bipartite if its vertices can be partitioned into
 * two sets, V1 and V2, so that connections only run between V1 and V2,
 * but not within V1 or within V2. The \p types parameter specifies
 * which vertex should be considered a member of one or the other
 * partition. The projection to V1 has vertex set V1, and two vertices
 * are connected if they have at least one common neighbour in V2.
 * The number of common neighbours is returned in \p multiplicity1,
 * if requested.
 *
 * \param graph The bipartite input graph. Directedness of the edges
 *   is ignored.
 * \param types Boolean vector giving the vertex types of the graph.
 * \param proj1 Pointer to an uninitialized graph object, the first
 *   projection will be created here. It a null pointer, then it is
 *   ignored, see also the \p probe1 argument.
 * \param proj2 Pointer to an uninitialized graph object, the second
 *   projection is created here, if it is not a null pointer. See also
 *   the \p probe1 argument.
 * \param multiplicity1 Pointer to a vector, or a null pointer. If not
 *   the latter, then the multiplicity of the edges is stored
 *   here. E.g. if there is an A-C-B and also an A-D-B triple in the
 *   bipartite graph (but no more X, such that A-X-B is also in the
 *   graph), then the multiplicity of the A-B edge in the projection
 *   will be 2.
 * \param multiplicity2 The same as \c multiplicity1, but for the
 *   other projection.
 * \param probe1 This argument can be used to specify the order of the
 *   projections in the resulting list. When it is non-negative, then
 *   it is considered as a vertex ID and the projection containing
 *   this vertex will be the first one in the result. Setting this
 *   argument to a non-negative value implies that \c proj1 must be
 *   a non-null pointer. If you don't care about the ordering of the
 *   projections, pass -1 here.
 * \return Error code.
 *
 * \sa \ref igraph_bipartite_projection_size() to calculate the number
 * of vertices and edges in the projections, without creating the
 * projection graphs themselves.
 *
 * Time complexity: O(|V|*d^2+|E|), |V| is the number of vertices, |E|
 * is the number of edges, d is the average (total) degree of the
 * graphs.
 */

igraph_error_t igraph_bipartite_projection(const igraph_t *graph,
                                const igraph_vector_bool_t *types,
                                igraph_t *proj1,
                                igraph_t *proj2,
                                igraph_vector_int_t *multiplicity1,
                                igraph_vector_int_t *multiplicity2,
                                igraph_int_t probe1) {

    const igraph_int_t no_of_nodes = igraph_vcount(graph);

    /* t1 is -1 if proj1 is omitted, it is 0 if it belongs to type zero,
       it is 1 if it belongs to type one. The same for t2 */
    int t1, t2;

    if (igraph_vector_bool_size(types) != no_of_nodes) {
        IGRAPH_ERROR("Invalid bipartite type vector length.", IGRAPH_EINVAL);
    }

    if (probe1 >= no_of_nodes) {
        IGRAPH_ERROR("No such vertex to probe.", IGRAPH_EINVAL);
    }

    if (probe1 >= 0 && !proj1) {
        IGRAPH_ERROR("`probe1' given, but `proj1' is a null pointer.", IGRAPH_EINVAL);
    }

    if (probe1 >= 0) {
        t1 = VECTOR(*types)[probe1];
        if (proj2) {
            t2 = 1 - t1;
        } else {
            t2 = -1;
        }
    } else {
        t1 = proj1 ? 0 : -1;
        t2 = proj2 ? 1 : -1;
    }

    if (proj1) {
        IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj1, t1, multiplicity1));
        IGRAPH_FINALLY(igraph_destroy, proj1);
    }

    if (proj2) {
        IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj2, t2, multiplicity2));
    }

    if (proj1) {
        IGRAPH_FINALLY_CLEAN(1); /* proj1 ownership change */
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_full_bipartite
 * \brief Creates a complete bipartite graph.
 *
 * A bipartite network contains two kinds of vertices and connections
 * are only possible between two vertices of different kind. There are
 * many natural examples, e.g. movies and actors as vertices and a
 * movie is connected to all participating actors, etc.
 *
 * </para><para>
 * igraph does not have direct support for bipartite networks, at
 * least not at the C language level. In other words the \type igraph_t
 * structure does not contain information about the vertex types.
 * The C functions for bipartite networks usually have an additional
 * input argument to graph, called \p types, a boolean vector giving
 * the vertex types.
 *
 * </para><para>
 * Most functions creating bipartite networks are able to create this
 * extra vector, you just need to supply an initialized boolean vector
 * to them.
 *
 * \param graph Pointer to an uninitialized graph object, the graph will be
 *   created here.
 * \param types Pointer to a boolean vector. If not a null pointer,
 *   then the vertex types will be stored here.
 * \param n1 Integer, the number of vertices of the first kind.
 * \param n2 Integer, the number of vertices of the second kind.
 * \param directed Boolean, whether to create a directed graph.
 * \param mode A constant that gives the type of connections for
 *   directed graphs. If \c IGRAPH_OUT, then edges point from vertices
 *   of the first kind to vertices of the second kind; if \c
 *   IGRAPH_IN, then the opposite direction is realized; if \c
 *   IGRAPH_ALL, then mutual edges will be created.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa \ref igraph_full() for non-bipartite complete graphs,
 * \ref igraph_full_multipartite() for complete multipartite graphs.
 */

igraph_error_t igraph_full_bipartite(igraph_t *graph,
                          igraph_vector_bool_t *types,
                          igraph_int_t n1, igraph_int_t n2,
                          igraph_bool_t directed,
                          igraph_neimode_t mode) {

    igraph_int_t no_of_nodes, no_of_edges;
    igraph_vector_int_t edges;
    igraph_int_t ptr;

    if (n1 < 0 || n2 < 0) {
        IGRAPH_ERROR("Invalid number of vertices for bipartite graph.", IGRAPH_EINVAL);
    }

    IGRAPH_SAFE_ADD(n1, n2, &no_of_nodes);

    if (!directed) {
        IGRAPH_SAFE_MULT(n1, n2, &no_of_edges);
    } else if (mode == IGRAPH_OUT || mode == IGRAPH_IN) {
        IGRAPH_SAFE_MULT(n1, n2, &no_of_edges);
    } else { /* mode==IGRAPH_ALL */
        IGRAPH_SAFE_MULT(n1, n2, &no_of_edges);
        IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges);
    }

    /* To ensure the size of the edges vector will not overflow. */
    if (no_of_edges > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Overflow in number of edges.", IGRAPH_EOVERFLOW);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);

    ptr = 0;

    if (!directed || mode == IGRAPH_OUT) {

        for (igraph_int_t i = 0; i < n1; i++) {
            for (igraph_int_t j = 0; j < n2; j++) {
                VECTOR(edges)[ptr++] = i;
                VECTOR(edges)[ptr++] = n1 + j;
            }
        }

    } else if (mode == IGRAPH_IN) {

        for (igraph_int_t i = 0; i < n1; i++) {
            for (igraph_int_t j = 0; j < n2; j++) {
                VECTOR(edges)[ptr++] = n1 + j;
                VECTOR(edges)[ptr++] = i;
            }
        }

    } else {

        for (igraph_int_t i = 0; i < n1; i++) {
            for (igraph_int_t j = 0; j < n2; j++) {
                VECTOR(edges)[ptr++] = i;
                VECTOR(edges)[ptr++] = n1 + j;
                VECTOR(edges)[ptr++] = n1 + j;
                VECTOR(edges)[ptr++] = i;
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_destroy, graph);

    if (types) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
        igraph_vector_bool_null(types);
        for (igraph_int_t i = n1; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = true;
        }
    }

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_create_bipartite
 * \brief Create a bipartite graph.
 *
 * This is a simple wrapper function to create a bipartite graph. It
 * does a little more than \ref igraph_create(), e.g. it checks that
 * the graph is indeed bipartite with respect to the given \p types
 * vector. If there is an edge connecting two vertices of the same
 * kind, then an error is reported.
 *
 * \param graph Pointer to an uninitialized graph object, the result is
 *   created here.
 * \param types Boolean vector giving the vertex types. The length of
 *   the vector defines the number of vertices in the graph.
 * \param edges Vector giving the edges of the graph. The highest
 *   vertex ID in this vector must be smaller than the length of the
 *   \p types vector.
 * \param directed Boolean, whether to create a directed graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \example examples/simple/igraph_bipartite_create.c
 */

igraph_error_t igraph_create_bipartite(igraph_t *graph, const igraph_vector_bool_t *types,
                            const igraph_vector_int_t *edges,
                            igraph_bool_t directed) {

    igraph_int_t no_of_nodes = igraph_vector_bool_size(types);
    igraph_int_t no_of_edges = igraph_vector_int_size(edges);
    igraph_int_t i;

    if (no_of_edges % 2 != 0) {
        IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVAL);
    }
    no_of_edges /= 2;

    if (! igraph_vector_int_isininterval(edges, 0, no_of_nodes-1)) {
        IGRAPH_ERROR("Invalid (negative or too large) vertex ID", IGRAPH_EINVVID);
    }

    /* Check bipartiteness */
    for (i = 0; i < no_of_edges * 2; i += 2) {
        igraph_int_t from = VECTOR(*edges)[i];
        igraph_int_t to = VECTOR(*edges)[i + 1];
        igraph_bool_t t1 = VECTOR(*types)[from];
        igraph_bool_t t2 = VECTOR(*types)[to];
        if ( (t1 && t2) || (!t1 && !t2) ) {
            IGRAPH_ERROR("Invalid edges, not a bipartite graph", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_empty(graph, no_of_nodes, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_edges(graph, edges, 0));

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_biadjacency
 * \brief Creates a bipartite graph from a bipartite adjacency matrix.
 *
 * A bipartite (or two-mode) graph contains two types of vertices and
 * edges always connect vertices of different types. A bipartite adjacency
 * matrix is an \em n x \em m matrix, \em n and \em m are the number of vertices
 * of the two types, respectively. Nonzero elements in the matrix denote
 * edges between the two corresponding vertices.
 *
 * </para><para>
 * This function can operate in two modes, depending on the
 * \p multiple argument. If it is \c false, then a single edge is
 * created for every non-zero element in the bipartite adjacency matrix. If
 * \p multiple is \c true, then as many edges are created between two
 * vertices as the corresponding matrix element. When \p multiple
 * is set to \c true, matrix elements should be whole numbers.
 * Otherwise their fractional part will be discarded.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param types Pointer to an initialized boolean vector, or a null
 *   pointer. If not a null pointer, then the vertex types are stored
 *   here. It is resized as needed.
 * \param biadjmatrix The bipartite adjacency matrix that serves as an input
 *   to this function.
 * \param directed Specifies whether to create an undirected or a directed
 *   graph.
 * \param mode Specifies the direction of the edges in a directed
 *   graph. If \c IGRAPH_OUT, then edges point from vertices
 *   of the first kind (corresponding to rows) to vertices of the
 *   second kind (corresponding to columns); if \c IGRAPH_IN,
 *   then the opposite direction is realized; if \c IGRAPH_ALL,
 *   then mutual edges will be created.
 * \param multiple Whether to interpret matrix entries as edge multiplicities,
 *   see details above.
 * \return Error code.
 *
 * Time complexity: O(n*m), the size of the bipartite adjacency matrix.
 */

igraph_error_t igraph_biadjacency(
        igraph_t *graph,
        igraph_vector_bool_t *types,
        const igraph_matrix_t *biadjmatrix,
        igraph_bool_t directed,
        igraph_neimode_t mode,
        igraph_bool_t multiple) {

    const igraph_int_t n1 = igraph_matrix_nrow(biadjmatrix);
    const igraph_int_t n2 = igraph_matrix_ncol(biadjmatrix);
    const igraph_int_t no_of_nodes = n1 + n2;
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (multiple) {

        if (n1 > 0 && n2 > 0 && igraph_matrix_min(biadjmatrix) < 0) {
            IGRAPH_ERRORF(
                "Bipartite adjacency matrix elements should be non-negative, found %g.",
                IGRAPH_EINVAL, igraph_matrix_min(biadjmatrix)
            );
        }

        for (igraph_int_t j = 0; j < n2; j++) {
            for (igraph_int_t i = 0; i < n1; i++) {
                igraph_int_t elem = MATRIX(*biadjmatrix, i, j);
                igraph_int_t from, to;

                if (elem == 0) {
                    continue;
                }

                if (mode == IGRAPH_IN) {
                    from = n1 + j;
                    to = i;
                } else {
                    from = i;
                    to = n1 + j;
                }

                if (mode != IGRAPH_ALL || !directed) {
                    for (igraph_int_t k = 0; k < elem; k++) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                    }
                } else {
                    for (igraph_int_t k = 0; k < elem; k++) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                    }
                }
            }
        }

    } else {

        for (igraph_int_t j = 0; j < n2; j++) {
            for (igraph_int_t i = 0; i < n1; i++) {
                igraph_int_t from, to;

                if (MATRIX(*biadjmatrix, i, j) != 0) {
                    if (mode == IGRAPH_IN) {
                        from = n1 + j;
                        to = i;
                    } else {
                        from = i;
                        to = n1 + j;
                    }
                    if (mode != IGRAPH_ALL || !directed) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                    } else {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                    }
                }
            }
        }

    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, graph);

    if (types) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
        igraph_vector_bool_null(types);
        for (igraph_int_t i = n1; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = true;
        }
    }

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_weighted_biadjacency
 * \brief Creates a bipartite graph from a weighted bipartite adjacency matrix.
 *
 * A bipartite (or two-mode) graph contains two types of vertices and
 * edges always connect vertices of different types. A bipartite adjacency
 * matrix is an \em n x \em m matrix, \em n and \em m are the number of vertices
 * of the two types, respectively. Nonzero elements in the matrix denote
 * edges between the two corresponding vertices.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param types Pointer to an initialized boolean vector, or a null
 *   pointer. If not a null pointer, then the vertex types are stored
 *   here. It is resized as needed.
 * \param weights Pointer to an initialized vector, the weights will be stored here.
 * \param biadjmatrix The bipartite adjacency matrix that serves as an input
 *   to this function.
 * \param directed Specifies whether to create an undirected or a directed
 *   graph.
 * \param mode Specifies the direction of the edges in a directed
 *   graph. If \c IGRAPH_OUT, then edges point from vertices
 *   of the first kind (corresponding to rows) to vertices of the
 *   second kind (corresponding to columns); if \c IGRAPH_IN,
 *   then the opposite direction is realized; if \c IGRAPH_ALL,
 *   then mutual edges will be created.
 * \return Error code.
 *
 * Time complexity: O(n*m), the size of the bipartite adjacency matrix.
 */

igraph_error_t igraph_weighted_biadjacency(
        igraph_t *graph,
        igraph_vector_bool_t *types,
        igraph_vector_t *weights,
        const igraph_matrix_t *biadjmatrix,
        igraph_bool_t directed,
        igraph_neimode_t mode) {

    const igraph_int_t n1 = igraph_matrix_nrow(biadjmatrix);
    const igraph_int_t n2 = igraph_matrix_ncol(biadjmatrix);
    const igraph_int_t no_of_nodes = n1 + n2;
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    igraph_vector_clear(weights);

    for (igraph_int_t j = 0; j < n2; j++) {
        for (igraph_int_t i = 0; i < n1; i++) {
            igraph_real_t weight = MATRIX(*biadjmatrix, i, j);
            igraph_int_t from, to;

            if (weight != 0) {
                if (mode == IGRAPH_IN) {
                    from = n1 + j;
                    to = i;
                } else {
                    from = i;
                    to = n1 + j;
                }
                if (mode != IGRAPH_ALL || !directed) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                    IGRAPH_CHECK(igraph_vector_push_back(weights, weight));
                } else {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                    IGRAPH_CHECK(igraph_vector_push_back(weights, weight));

                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                    IGRAPH_CHECK(igraph_vector_push_back(weights, weight));
                }
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, graph);

    if (types) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
        igraph_vector_bool_null(types);
        for (igraph_int_t i = n1; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = true;
        }
    }

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_biadjacency
 * \brief Converts a bipartite graph into a bipartite adjacency matrix.
 *
 * In a bipartite adjacency matrix \c A, element <code>A_ij</code>
 * gives the number of edges between the <code>i</code>th vertex of the
 * first partition and the <code>j</code>th vertex of the second partition.
 *
 * </para><para>
 * If the graph contains edges within the same partition, this function
 * issues a warning.
 *
 * \param graph The input graph, edge directions are ignored.
 * \param types Boolean vector containing the vertex types. Vertices belonging
 *   to the first partition have type \c false, the one in the second
 *   partition type \c true.
 * \param weights A vector specifying a weight for each edge or \c NULL.
 *   If \c NULL, all edges are assumed to have weight 1.
 * \param res Pointer to an initialized matrix, the result is stored
 *   here. An element of the matrix gives the number of edges
 *   (irrespectively of their direction), or sum of edge weights,
 *   between the two corresponding vertices. The rows will correspond
 *   to vertices with type \c false, the columns correspond to vertices
 *   with type \c true.
 * \param row_ids Pointer to an initialized vector or \c NULL.
 *   If not a null pointer, then the IDs of vertices with type \c false
 *   are stored here, with the same ordering as the rows of the
 *   biadjacency matrix.
 * \param col_ids Pointer to an initialized vector or \c NULL.
 *   If not a null pointer, then the IDs of vertices with type \c true
 *   are stored here, with the same ordering as the columns of the
 *   biadjacency matrix.
 * \return Error code.
 *
 * Time complexity: O(|E|) where |E| is the number of edges.
 *
 * \sa \ref igraph_biadjacency() for the opposite operation.
 */

igraph_error_t igraph_get_biadjacency(
    const igraph_t *graph, const igraph_vector_bool_t *types,
    const igraph_vector_t *weights,
    igraph_matrix_t *res, igraph_vector_int_t *row_ids,
    igraph_vector_int_t *col_ids
) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_int_t n1 = 0, n2 = 0;
    igraph_int_t ignored_edges = 0;
    igraph_vector_int_t perm;

    if (igraph_vector_bool_size(types) != no_of_nodes) {
        IGRAPH_ERRORF("Vertex type vector size (%" IGRAPH_PRId ") not equal to number of vertices (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_bool_size(types), no_of_nodes);
    }

    if (weights) {
        if (igraph_vector_size(weights) != no_of_edges) {
            IGRAPH_ERRORF("Edge weight vector size (%" IGRAPH_PRId ") not equal to number of edges (%" IGRAPH_PRId ").",
                          IGRAPH_EINVAL, igraph_vector_size(weights), no_of_edges);
        }
    }

    for (igraph_int_t i = 0; i < no_of_nodes; i++) {
        n1 += VECTOR(*types)[i] == false ? 1 : 0;
    }
    n2 = no_of_nodes - n1;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&perm, no_of_nodes);

    for (igraph_int_t i = 0, p1 = 0, p2 = n1; i < no_of_nodes; i++) {
        VECTOR(perm)[i] = VECTOR(*types)[i] ? p2++ : p1++;
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, n1, n2));
    igraph_matrix_null(res);
    for (igraph_int_t i = 0; i < no_of_edges; i++) {
        igraph_int_t from = IGRAPH_FROM(graph, i);
        igraph_int_t to = IGRAPH_TO(graph, i);
        igraph_int_t from2 = VECTOR(perm)[from];
        igraph_int_t to2 = VECTOR(perm)[to];
        if (VECTOR(*types)[from] == VECTOR(*types)[to]) {
            ignored_edges++;
        } else if (! VECTOR(*types)[from]) {
            MATRIX(*res, from2, to2 - n1) += weights ? VECTOR(*weights)[i] : 1;
        } else {
            MATRIX(*res, to2, from2 - n1) += weights ? VECTOR(*weights)[i] : 1;
        }
    }

    if (ignored_edges > 0) {
        IGRAPH_WARNINGF("%" IGRAPH_PRId " edges running within partitions were ignored.", ignored_edges);
    }

    if (row_ids) {
        IGRAPH_CHECK(igraph_vector_int_resize(row_ids, n1));
    }
    if (col_ids) {
        IGRAPH_CHECK(igraph_vector_int_resize(col_ids, n2));
    }
    if (row_ids || col_ids) {
        for (igraph_int_t i = 0; i < no_of_nodes; i++) {
            if (! VECTOR(*types)[i]) {
                if (row_ids) {
                    igraph_int_t i2 = VECTOR(perm)[i];
                    VECTOR(*row_ids)[i2] = i;
                }
            } else {
                if (col_ids) {
                    igraph_int_t i2 = VECTOR(perm)[i];
                    VECTOR(*col_ids)[i2 - n1] = i;
                }
            }
        }
    }

    igraph_vector_int_destroy(&perm);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_bipartite
 * \brief Check whether a graph is bipartite.
 *
 * This function checks whether a graph is bipartite. It tries
 * to find a mapping that gives a possible division of the vertices into two
 * classes, such that no two vertices of the same class are connected by an
 * edge.
 *
 * </para><para>
 * The existence of such a mapping is equivalent of having no circuits of
 * odd length in the graph. A graph with loop edges cannot be bipartite.
 *
 * </para><para>
 * Note that the mapping is not necessarily unique, e.g. if the graph has
 * at least two components, then the vertices in the separate components
 * can be mapped independently.
 *
 * \param graph The input graph.
 * \param res Pointer to a boolean, the result is stored here.
 * \param types Pointer to an initialized boolean vector, or a null
 *   pointer. If not a null pointer and a mapping was found, then it
 *   is stored here. If not a null pointer, but no mapping was found,
 *   the contents of this vector is invalid.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \sa igraph_is_bipartite_coloring() to determine if all edges connect
 * vertices of different types, given a specific type vector.
 */

igraph_error_t igraph_is_bipartite(const igraph_t *graph,
                        igraph_bool_t *res,
                        igraph_vector_bool_t *types) {

    /* We basically do a breadth first search and label the
       vertices along the way. We stop as soon as we can find a
       contradiction.

       In the 'seen' vector 0 means 'not seen yet', 1 means type 1,
       2 means type 2.
    */

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_char_t seen;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t neis;
    igraph_bool_t bi = true;

    /* Shortcut: Graphs with self-loops are not bipartite. */
    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_LOOP) &&
        igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_LOOP)) {
        if (res) {
            *res = false;
        }
        return IGRAPH_SUCCESS;
    }

    /* Shortcut: If the type vector is not requested, and the graph is a forest
     * we can immediately return with the result that the graph is bipartite. */
    if (! types &&
        igraph_i_property_cache_has(graph, IGRAPH_PROP_IS_FOREST) &&
        igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_IS_FOREST)) {
        if (res) {
            *res = true;
        }
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&seen, no_of_nodes);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    for (igraph_int_t i = 0; bi && i < no_of_nodes; i++) {

        if (VECTOR(seen)[i]) {
            continue;
        }

        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, i));
        VECTOR(seen)[i] = 1;

        while (bi && !igraph_dqueue_int_empty(&Q)) {
            igraph_int_t n, j;
            igraph_int_t actnode = igraph_dqueue_int_pop(&Q);
            char acttype = VECTOR(seen)[actnode];

            IGRAPH_CHECK(igraph_neighbors(
                graph, &neis, actnode, IGRAPH_ALL, IGRAPH_LOOPS, IGRAPH_MULTIPLE
            ));
            n = igraph_vector_int_size(&neis);
            for (j = 0; j < n; j++) {
                igraph_int_t nei = VECTOR(neis)[j];
                if (VECTOR(seen)[nei]) {
                    char neitype = VECTOR(seen)[nei];
                    if (neitype == acttype) {
                        bi = false;
                        break;
                    }
                } else {
                    VECTOR(seen)[nei] = 3 - acttype;
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, nei));
                }
            }
        }
    }

    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&Q);
    IGRAPH_FINALLY_CLEAN(2);

    /* Set the cache: A graph that is not bipartite has
     * an odd-length cycle, therefore it cannot be a forest. */
    if (! bi) {
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_IS_FOREST, false);
    }

    if (res) {
        *res = bi;
    }

    if (types && bi) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
        for (igraph_int_t i = 0; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = VECTOR(seen)[i] - 1;
        }
    }

    igraph_vector_char_destroy(&seen);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


static igraph_error_t bipartite_iea_game(
        igraph_t *graph,
        igraph_int_t n1, igraph_int_t n2,
        igraph_int_t m,
        igraph_bool_t directed, igraph_neimode_t mode) {

    igraph_vector_int_t edges;
    igraph_int_t n = n1 + n2; /* overflow checked by caller */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, m * 2));

    for (igraph_int_t i = 0; i < m; i++) {
        igraph_int_t to, from;

        to = RNG_INTEGER(n1, n - 1);
        from = RNG_INTEGER(0, n1 - 1);

        /* flip unconditionally for IGRAPH_IN,
         * or with probability 0.5 for IGRAPH_ALL */
        if (mode == IGRAPH_IN || (mode == IGRAPH_ALL && RNG_BOOL())) {
            igraph_vector_int_push_back(&edges, to);
            igraph_vector_int_push_back(&edges, from);
        } else {
            igraph_vector_int_push_back(&edges, from);
            igraph_vector_int_push_back(&edges, to);
        }

    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t bipartite_gnm_multi(
        igraph_t *graph,
        igraph_int_t n1, igraph_int_t n2,
        igraph_int_t m,
        igraph_bool_t directed, igraph_neimode_t mode) {

    /* See igraph_erdos_renyi_game_gnm() for how the sampling works. */

    igraph_vector_int_t edges;
    igraph_int_t nrow, ncol;
    igraph_int_t last;
    igraph_int_t offset1 = 0, offset2 = n1;
    igraph_int_t n = n1 + n2; /* overflow checked by caller */

    /* The larger partition is associated with columns, the smaller
     * with rows. This setup helps avoid integer overflow. We swap
     * n1 and n2 so that n1 is smaller. */
    if (n1 > n2) {
        igraph_int_t tmp = n1;
        n1 = n2;
        n2 = tmp;

        offset1 = n2; offset2 = 0;

        mode = IGRAPH_REVERSE_MODE(mode);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2*m);

    if (!directed || mode != IGRAPH_ALL) {
        nrow = n1;
        ncol = n2;
        last = ncol-1;
        for (igraph_int_t i=0; i < m; i++) {
            while (true) {
                igraph_int_t r = RNG_INTEGER(0, nrow-1);
                igraph_int_t c = RNG_INTEGER(0, ncol-1);

                if (r >= n1) {
                    igraph_int_t j = (r - n1) * ncol + c;
                    if (j >= i) continue; /* rejection sampling */
                    VECTOR(edges)[2*i]   = VECTOR(edges)[2*j];
                    VECTOR(edges)[2*i+1] = VECTOR(edges)[2*j+1];
                } else {
                    if (directed && mode == IGRAPH_IN) {
                        VECTOR(edges)[2*i]   = c + offset2;
                        VECTOR(edges)[2*i+1] = r + offset1;
                    } else {
                        VECTOR(edges)[2*i]   = r + offset1;
                        VECTOR(edges)[2*i+1] = c + offset2;
                    }
                }

                last += 1;
                if (last >= ncol) {
                    last -= ncol;
                    nrow++;
                }

                break;
            }
        }
    } else /* directed, mutual allowed */ {
        nrow = 2*n1;
        ncol = n2;
        last = ncol-1;
        for (igraph_int_t i=0; i < m; i++) {
            while (true) {
                igraph_int_t r = RNG_INTEGER(0, nrow-1);
                igraph_int_t c = RNG_INTEGER(0, ncol-1);

                if (r >= 2*n1) {
                    igraph_int_t j = (r - 2*n1) * ncol + c;
                    if (j >= i) continue; /* rejection sampling */
                    VECTOR(edges)[2*i]   = VECTOR(edges)[2*j];
                    VECTOR(edges)[2*i+1] = VECTOR(edges)[2*j+1];
                } else {
                    if (r < n1) {
                        VECTOR(edges)[2*i]   = r + offset1;
                        VECTOR(edges)[2*i+1] = c + offset2;
                    } else {
                        VECTOR(edges)[2*i]   = c + offset2;
                        VECTOR(edges)[2*i+1] = r - n1 + offset1;
                    }
                }

                last += 1;
                if (last >= ncol) {
                    last -= ncol;
                    nrow++;
                }

                break;
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t bipartite_gnm_simple(
        igraph_t *graph,
        igraph_int_t n1, igraph_int_t n2,
        igraph_int_t m,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_bool_t edge_labeled) {

    igraph_vector_int_t edges;
    igraph_vector_t s;
    igraph_real_t n1_real = (igraph_real_t) n1, n2_real = (igraph_real_t) n2; /* for floating-point operations */
    igraph_int_t n = n1 + n2; /* overflow checked by caller */
    igraph_real_t maxedges;
    int iter = 0;

    if (!directed || mode != IGRAPH_ALL) {
        maxedges = n1_real * n2_real;
    } else {
        maxedges = 2.0 * n1_real * n2_real;
    }

    if (m > maxedges) {
        IGRAPH_ERROR("Too many edges requested compared to the number of vertices.", IGRAPH_EINVAL);
    }

    if (maxedges == m && ! edge_labeled) {
        /* TODO: Cannot use igraph_full_bipartite() when edge_labeled as we must shuffle edges. */
        IGRAPH_CHECK(igraph_full_bipartite(graph, NULL, n1, n2, directed, mode));
    } else {
        igraph_int_t to, from;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_i_random_sample_real(&s, 0, maxedges - 1, m));
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, m * 2));

        for (igraph_int_t i = 0; i < m; i++) {
            if (!directed || mode != IGRAPH_ALL) {
                to = trunc(VECTOR(s)[i] / n1_real);
                from = VECTOR(s)[i] - to * n1_real;
                to += n1;
            } else {
                igraph_real_t n1n2 = n1_real * n2_real;
                if (VECTOR(s)[i] < n1n2) {
                    to = trunc(VECTOR(s)[i] / n1_real);
                    from = VECTOR(s)[i] - to * n1_real;
                    to += n1;
                } else {
                    to = trunc((VECTOR(s)[i] - n1n2) / n2_real);
                    from = VECTOR(s)[i] - n1n2 - to * n2_real;
                    from += n1;
                }
            }

            if (mode != IGRAPH_IN) {
                igraph_vector_int_push_back(&edges, from); /* reserved */
                igraph_vector_int_push_back(&edges, to); /* reserved */
            } else {
                igraph_vector_int_push_back(&edges, to); /* reserved */
                igraph_vector_int_push_back(&edges, from); /* reserved */
            }

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }

        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);

        if (edge_labeled) {
            IGRAPH_CHECK(igraph_i_vector_int_shuffle_pairs(&edges));
        }

        IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));

        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bipartite_game_gnm
 * \brief Generate a random bipartite graph with a fixed number of edges.
 *
 * The <code>G(n1, n2, m)</code> model uniformly samples bipartite graphs with
 * \p n1 bottom vertices and \p n2 top vertices, and precisely \p m edges.
 *
 * \param graph Pointer to an uninitialized igraph graph, the result
 *    is stored here.
 * \param types Pointer to an initialized boolean vector, or a null
 *    pointer. If not a null pointer, then the vertex types are stored
 *    here. Bottom vertices come first, \p n1 of them, then \p n2 top
 *    vertices.
 * \param n1 The number of bottom vertices.
 * \param n2 The number of top vertices.
 * \param m The number of edges.
 * \param directed Boolean, whether to generate a directed graph. See
 *     also the \p mode argument.
 * \param mode Specifies how to direct the edges in directed
 *     graphs. If it is \c IGRAPH_OUT, then directed edges point from
 *     bottom vertices to top vertices. If it is \c IGRAPH_IN, edges
 *     point from top vertices to bottom vertices. \c IGRAPH_OUT and
 *     \c IGRAPH_IN do not generate mutual edges. If this argument is
 *     \c IGRAPH_ALL, then each edge direction is considered
 *     independently and mutual edges might be generated. This
 *     argument is ignored for undirected graphs.
* \param allowed_edge_types The types of edges to allow in the graph.
 *        \clist
 *          \cli IGRAPH_SIMPLE_SW
 *          simple graph (i.e. no multi-edges allowed).
 *          \cli IGRAPH_MULTI_SW
 *          multi-edges are allowed
 *        \endclist
 * \param edge_labeled If true, the sampling is done uniformly from the set
 *     of ordered edge lists. See \ref igraph_bipartite_iea_game() for more
 *     information. Set this to \c false to select the classic Erdős-Rényi model.
 *     The constants \c IGRAPH_EDGE_UNLABELED and \c IGRAPH_EDGE_LABELED
 *     may be used instead of \c false and \c true for better readability.
 * \return Error code.
 *
 * \sa \ref igraph_erdos_renyi_game_gnm() for the unipartite version,
 * \ref igraph_bipartite_game_gnp() for the <code>G(n1, n2, p)</code>
 * model.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

igraph_error_t igraph_bipartite_game_gnm(
        igraph_t *graph,
        igraph_vector_bool_t *types,
        igraph_int_t n1, igraph_int_t n2, igraph_int_t m,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_bool_t edge_labeled) {

    igraph_int_t n;
    igraph_bool_t loops, multiple;

    if (n1 < 0 || n2 < 0) {
        IGRAPH_ERROR("Invalid number of vertices for bipartite G(n,m) model.", IGRAPH_EINVAL);
    }
    if (m < 0 || m > IGRAPH_ECOUNT_MAX) {
        IGRAPH_ERROR("Invalid number of edges for bipartite G(n,m) model.", IGRAPH_EINVAL);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for bipartite G(n,m) model.", IGRAPH_EINVAL);
    }

    /* Bipartite graphs cannot have self-loops. We ignore them. */
    IGRAPH_CHECK(igraph_i_edge_type_to_loops_multiple(allowed_edge_types, &loops, &multiple));

    IGRAPH_SAFE_ADD(n1, n2, &n); /* overflow check */

    if (types) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, n));
        igraph_vector_bool_null(types);
        for (igraph_int_t i = n1; i < n; i++) {
            VECTOR(*types)[i] = true;
        }
    }

    if (m == 0 || n1 == 0 || n2 == 0) {
        if (m > 0) {
            IGRAPH_ERROR("Too many edges requested compared to the number of vertices.", IGRAPH_EINVAL);
        }
        return igraph_empty(graph, n, directed);
    } else if (multiple) {
        if (edge_labeled) {
            return bipartite_iea_game(graph, n1, n2, m, directed, mode);
        } else {
            return bipartite_gnm_multi(graph, n1, n2, m, directed, mode);
        }
    } else {
        return bipartite_gnm_simple(graph, n1, n2, m, directed, mode, edge_labeled);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bipartite_iea_game
 * \brief Generates a random bipartite multigraph through independent edge assignment.
 *
 * \experimental
 *
 * This model generates random multigraphs with \p n1 bottom vertices,
 * \p n2 top vertices and \p m edges through independent edge assignment (IEA).
 * Each of the \p m edges is assigned uniformly at random to a vertex pair,
 * independently of each other.
 *
 * </para><para>
 * This model does not sample multigraphs uniformly. Undirected graphs are
 * generated with probability proportional to
 *
 * </para><para>
 * <code>(prod_(i&lt;j) A_ij !)^(-1)</code>,
 *
 * </para><para>
 * where \c A denotes the adjacency matrix. The corresponding  expression for
 * directed graphs is
 *
 * </para><para>
 * <code>(prod_(i,j) A_ij !)^(-1)</code>.
 *
 * </para><para>
 * Thus the probability of all simple graphs (which only have 0s and 1s in the
 * adjacency matrix) is the same, while that of non-simple ones depends on
 * their edge and self-loop multiplicities.
 *
 * \param graph Pointer to an uninitialized igraph graph, the result
 *    is stored here.
 * \param types Pointer to an initialized boolean vector, or a \c NULL
 *    pointer. If not \c NULL, then the vertex types are stored
 *    here. Bottom vertices come first, \p n1 of them, then \p n2 top
 *    vertices.
 * \param n1 The number of bottom vertices.
 * \param n2 The number of top vertices.
 * \param m The number of edges.
 * \param directed Whether to generate a directed graph. See
 *     also the \p mode argument.
 * \param mode Specifies how to direct the edges in directed
 *     graphs. If it is \c IGRAPH_OUT, then directed edges point from
 *     bottom vertices to top vertices. If it is \c IGRAPH_IN, edges
 *     point from top vertices to bottom vertices. \c IGRAPH_OUT and
 *     \c IGRAPH_IN do not generate mutual edges. If this argument is
 *     \c IGRAPH_ALL, then each edge direction is considered
 *     independently and mutual edges might be generated. This
 *     argument is ignored for undirected graphs.
 * \return Error code.
 *
 * \sa \ref igraph_iea_game() for the unipartite version;
 * \ref igraph_bipartite_game_gnm() to uniformly sample bipartite graphs
 * with a given number of vertices and edges.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

igraph_error_t igraph_bipartite_iea_game(
        igraph_t *graph, igraph_vector_bool_t *types,
        igraph_int_t n1, igraph_int_t n2, igraph_int_t m,
        igraph_bool_t directed, igraph_neimode_t mode) {

    return igraph_bipartite_game_gnm(
        graph, types, n1, n2, m, directed, mode,
        IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);
}


static igraph_error_t bipartite_gnp_edge_labeled(
        igraph_t *graph,
        igraph_int_t n1, igraph_int_t n2, igraph_real_t p,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_bool_t multiple) {

    if (multiple) {
        igraph_real_t maxedges;

        if (!directed || mode != IGRAPH_ALL) {
            maxedges = (igraph_real_t) n1 * (igraph_real_t) n2;
        } else {
            maxedges = 2.0 * (igraph_real_t) n1 * (igraph_real_t) n2;
        }

        igraph_real_t m;
        do {
            m = RNG_GEOM( 1.0 / (1.0 + maxedges * p) );
        } while (m > (igraph_real_t) IGRAPH_INTEGER_MAX);

        return bipartite_iea_game(graph, n1, n2, m, directed, mode);
    } else {
        IGRAPH_ERROR("The edge-labeled bipartite G(n,p) model is not yet implemented for graphs without multi-edges.",
                     IGRAPH_UNIMPLEMENTED);
    }
}

/* This implementation is used only with very large vertex counts, when the
 * default implementation would fail due to overflow. While this version
 * avoids overflow and uses less memory, it is also slower than the default
 * implementation.
 *
 * This function expects that when multiple=true, the p parameter has already
 * been transformed by p = p / (1 + p). This is currently done by the caller.
 */
static igraph_error_t gnp_bipartite_large(
        igraph_t *graph,
        igraph_int_t n1, igraph_int_t n2,
        igraph_real_t p,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_bool_t multiple,
        igraph_int_t ecount_estimate) {

    igraph_vector_int_t edges;
    int iter = 0;

    /* Necessitated by floating point arithmetic used in the implementation. */
    if (n1 >= IGRAPH_MAX_EXACT_REAL || n2 >= IGRAPH_MAX_EXACT_REAL) {
        IGRAPH_ERROR("Number of vertices is too large.", IGRAPH_EOVERFLOW);
    }

    if (ecount_estimate > IGRAPH_ECOUNT_MAX) {
        ecount_estimate = IGRAPH_ECOUNT_MAX;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, 2*ecount_estimate));

    for (igraph_int_t i = 0; i < n1; i++) {
        igraph_int_t j = 0;

        while (true) {
            igraph_real_t gap = RNG_GEOM(p);

            if (gap >= n2 - j) {
                break;
            }

            j += gap;

            if (!directed) {
                /* Undirected graph */
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j + n1));
            } else if (mode == IGRAPH_IN) {
                /* Incoming edges */
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j + n1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            } else if (mode == IGRAPH_OUT) {
                /* Outgoing edges */
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j + n1));
            } else {
                /* Both directions for IGRAPH_ALL */
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j + n1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j + n1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            }

            j += ! multiple; /* 1 for simple graph, 0 for multigraph */

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }
    }

    /* n1 + n2 has already been checked for overflow in the caller function. */
    IGRAPH_CHECK(igraph_create(graph, &edges, n1 + n2, directed));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bipartite_game_gnp
 * \brief Generates a random bipartite graph with a fixed connection probability.
 *
 * In the <code>G(n1, n2, p)</code> model, every possible edge between the \p n1
 * bottom vertices and \p n2 top vertices is realized independently with
 * probability \p p. This is equivalent to a maximum entropy model with
 * a constraint on the \em expected total edge count. This view allows
 * a multigraph extension, in which case \p is interpreted as the expected
 * number of edges between any vertex pair. See \ref igraph_erdos_renyi_game_gnp()
 * for more details.
 *
 * \param graph Pointer to an uninitialized igraph graph, the result
 *    is stored here.
 * \param types Pointer to an initialized boolean vector, or a null
 *    pointer. If not \c NULL, then the vertex types are stored
 *    here. Bottom vertices come first, \p n1 of them, then \p n2 top
 *    vertices.
 * \param n1 The number of bottom vertices.
 * \param n2 The number of top vertices.
 * \param p The expected number of edges between any vertex pair.
 *    When multi-edges are disallowed, this is equivalent to the probability
 *    of having a connection between any two vertices.
 * \param directed Whether to generate a directed graph. See also
 *     the \p mode argument.
 * \param mode Specifies how to direct the edges in directed
 *     graphs. If it is \c IGRAPH_OUT, then directed edges point from
 *     bottom vertices to top vertices. If it is \c IGRAPH_IN, edges
 *     point from top vertices to bottom vertices. \c IGRAPH_OUT and
 *     \c IGRAPH_IN do not generate mutual edges. If this argument is
 *     \c IGRAPH_ALL, then each edge direction is considered
 *     independently and mutual edges might be generated. This
 *     argument is ignored for undirected graphs.
* \param allowed_edge_types The types of edges to allow in the graph.
 *        \clist
 *          \cli IGRAPH_SIMPLE_SW
 *          simple graph (i.e. no multi-edges allowed).
 *          \cli IGRAPH_MULTI_SW
 *          multi-edges are allowed
 *        \endclist
 * \param edge_labeled If true, the model is defined over the set of ordered
 *     edge lists, i.e. over the set of edge-labeled graphs. Set it to
 *     \c false to select the classic bipartite Erdős-Rényi model.
 *     The constants \c IGRAPH_EDGE_UNLABELED and \c IGRAPH_EDGE_LABELED
 *     may be used instead of \c false and \c true for better readability.
 * \return Error code.
 *
 * \sa \ref igraph_erdos_renyi_game_gnp() for the unipartite version,
 * \ref igraph_bipartite_game_gnm() for the <code>G(n1, n2, m)</code> model.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

igraph_error_t igraph_bipartite_game_gnp(
        igraph_t *graph,
        igraph_vector_bool_t *types,
        igraph_int_t n1, igraph_int_t n2, igraph_real_t p,
        igraph_bool_t directed, igraph_neimode_t mode,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_bool_t edge_labeled) {

    igraph_vector_int_t edges;
    igraph_vector_t s;
    igraph_int_t n;
    igraph_real_t n1_real = (igraph_real_t) n1, n2_real = (igraph_real_t) n2; /* for floating-point operations */
    igraph_bool_t loops, multiple;
    int iter = 0;

    if (n1 < 0 || n2 < 0) {
        IGRAPH_ERROR("Invalid number of vertices for bipartite G(n,p) model.", IGRAPH_EINVAL);
    }

    /* Bipartite graphs cannot have self-loops. We ignore them. */
    IGRAPH_CHECK(igraph_i_edge_type_to_loops_multiple(allowed_edge_types, &loops, &multiple));

    if (multiple) {
        if (p < 0.0) {
            IGRAPH_ERROR(
                "Invalid expected edge multiplicity given for "
                "bipartite G(n,p) multigraph model.",
                IGRAPH_EINVAL);
        }
    } else {
        if (p < 0.0 || p > 1.0) {
            IGRAPH_ERROR(
                "Invalid connection probability given for bipartite G(n,p) model.",
                IGRAPH_EINVAL);
        }
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode for bipartite G(n,p) model.", IGRAPH_EINVAL);
    }

    IGRAPH_SAFE_ADD(n1, n2, &n);

    if (types) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, n));
        igraph_vector_bool_null(types);
        for (igraph_int_t i = n1; i < n; i++) {
            VECTOR(*types)[i] = true;
        }
    }

    if (edge_labeled) {
        return bipartite_gnp_edge_labeled(graph, n1, n2, p, directed, mode, multiple);
    }

    if (multiple) {
        /* Convert the expected edge count to the appropriate probability parameter
         * of the geometric distribution when sampling lengths of runs of 0s in the
         * adjacency matrix. */
        p = p / (1 + p);
    }

    if (p == 0 || n1 == 0 || n2 == 0) {
        IGRAPH_CHECK(igraph_empty(graph, n, directed));
    } else if (p == 1.0) {
        IGRAPH_CHECK(igraph_full_bipartite(graph, types, n1, n2, directed, mode));
    } else {

        igraph_int_t to, from, slen;
        igraph_real_t maxedges, last;
        igraph_int_t ecount_estimate;

        if (!directed || mode != IGRAPH_ALL) {
            maxedges = n1_real * n2_real;
        } else {
            maxedges = 2.0 * n1_real * n2_real;
        }

        IGRAPH_CHECK(igraph_i_safe_floor(maxedges * p * 1.1, &ecount_estimate));

        if (maxedges > IGRAPH_MAX_EXACT_REAL) {
            /* Use a slightly slower, but overflow-free implementation. */
            return gnp_bipartite_large(graph, n1, n2, p, directed, mode, multiple, ecount_estimate);
        }

        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, ecount_estimate));

        last = RNG_GEOM(p);
        while (last < maxedges) {
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(p);
            last += ! multiple; /* 1 for simple graph, 0 for multigraph */
            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }

        slen = igraph_vector_size(&s);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, slen * 2));

        for (igraph_int_t i = 0; i < slen; i++) {
            if (!directed || mode != IGRAPH_ALL) {
                to = trunc(VECTOR(s)[i] / n1_real);
                from = VECTOR(s)[i] - to * n1_real;
                to += n1;
            } else {
                igraph_real_t n1n2 = n1_real * n2_real;
                if (VECTOR(s)[i] < n1n2) {
                    to = trunc(VECTOR(s)[i] / n1_real);
                    from = VECTOR(s)[i] - to * n1_real;
                    to += n1;
                } else {
                    to = trunc((VECTOR(s)[i] - n1n2) / n2_real);
                    from = VECTOR(s)[i] - n1n2 - to * n2_real;
                    from += n1;
                }
            }

            if (mode != IGRAPH_IN) {
                igraph_vector_int_push_back(&edges, from); /* reserved */
                igraph_vector_int_push_back(&edges, to); /* reserved */
            } else {
                igraph_vector_int_push_back(&edges, to); /* reserved */
                igraph_vector_int_push_back(&edges, from); /* reserved */
            }

            IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 14);
        }

        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);

        IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
