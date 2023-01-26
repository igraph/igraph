/* -*- mode: C -*-  */
/*
   IGraph library.
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
#include "igraph_nongraph.h"

#include "graph/attributes.h"
#include "math/safe_intop.h"

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
 * \param vcount1 Pointer to an \c igraph_integer_t, the number of
 *     vertices in the first projection is stored here. May be \c NULL
 *     if not needed.
 * \param ecount1 Pointer to an \c igraph_integer_t, the number of
 *     edges in the first projection is stored here. May be \c NULL
 *     if not needed.
 * \param vcount2 Pointer to an \c igraph_integer_t, the number of
 *     vertices in the second projection is stored here. May be \c NULL
 *     if not needed.
 * \param ecount2 Pointer to an \c igraph_integer_t, the number of
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
                                     igraph_integer_t *vcount1,
                                     igraph_integer_t *ecount1,
                                     igraph_integer_t *vcount2,
                                     igraph_integer_t *ecount2) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t vc1 = 0, ec1 = 0, vc2 = 0, ec2 = 0;
    igraph_adjlist_t adjlist;
    igraph_vector_int_t added;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&added, no_of_nodes);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    for (i = 0; i < no_of_nodes; i++) {
        igraph_vector_int_t *neis1;
        igraph_integer_t neilen1, j;
        igraph_integer_t *ecptr;
        if (VECTOR(*types)[i]) {
            vc2++;
            ecptr = &ec2;
        } else {
            vc1++;
            ecptr = &ec1;
        }
        neis1 = igraph_adjlist_get(&adjlist, i);
        neilen1 = igraph_vector_int_size(neis1);
        for (j = 0; j < neilen1; j++) {
            igraph_integer_t k, neilen2, nei = VECTOR(*neis1)[j];
            igraph_vector_int_t *neis2 = igraph_adjlist_get(&adjlist, nei);
            if (IGRAPH_UNLIKELY(VECTOR(*types)[i] == VECTOR(*types)[nei])) {
                IGRAPH_ERROR("Non-bipartite edge found in bipartite projection",
                             IGRAPH_EINVAL);
            }
            neilen2 = igraph_vector_int_size(neis2);
            for (k = 0; k < neilen2; k++) {
                igraph_integer_t nei2 = VECTOR(*neis2)[k];
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

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t i, j, k;
    igraph_integer_t remaining_nodes = 0;
    igraph_vector_int_t vertex_perm, vertex_index;
    igraph_vector_int_t edges;
    igraph_adjlist_t adjlist;
    igraph_vector_int_t *neis1, *neis2;
    igraph_integer_t neilen1, neilen2;
    igraph_vector_int_t added;
    igraph_vector_t mult;

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
    IGRAPH_VECTOR_INIT_FINALLY(&mult, multiplicity ? no_of_nodes : 1);
    if (multiplicity) {
        igraph_vector_int_clear(multiplicity);
    }

    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] == which) {
            VECTOR(vertex_index)[i] = ++remaining_nodes;
            igraph_vector_int_push_back(&vertex_perm, i);
        }
    }

    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*types)[i] == which) {
            igraph_integer_t new_i = VECTOR(vertex_index)[i] - 1;
            igraph_integer_t iedges = 0;
            neis1 = igraph_adjlist_get(&adjlist, i);
            neilen1 = igraph_vector_int_size(neis1);
            for (j = 0; j < neilen1; j++) {
                igraph_integer_t nei = VECTOR(*neis1)[j];
                if (IGRAPH_UNLIKELY(VECTOR(*types)[i] == VECTOR(*types)[nei])) {
                    IGRAPH_ERROR("Non-bipartite edge found in bipartite projection",
                                 IGRAPH_EINVAL);
                }
                neis2 = igraph_adjlist_get(&adjlist, nei);
                neilen2 = igraph_vector_int_size(neis2);
                for (k = 0; k < neilen2; k++) {
                    igraph_integer_t nei2 = VECTOR(*neis2)[k], new_nei2;
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
                        new_nei2 = VECTOR(vertex_index)[nei2] - 1;
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, new_nei2));
                    }
                }
            }
            if (multiplicity) {
                /* OK, we need to go through all the edges added for vertex new_i
                   and check their multiplicity */
                igraph_integer_t now = igraph_vector_int_size(&edges);
                igraph_integer_t from = now - iedges * 2;
                for (j = from; j < now; j += 2) {
                    igraph_integer_t nei2 = VECTOR(edges)[j + 1];
                    igraph_integer_t new_nei2 = VECTOR(vertex_index)[nei2] - 1;
                    igraph_integer_t m = VECTOR(mult)[nei2];
                    VECTOR(edges)[j + 1] = new_nei2;
                    IGRAPH_CHECK(igraph_vector_int_push_back(multiplicity, m));
                }
            }
        } /* if VECTOR(*type)[i] == which */
    }

    igraph_vector_destroy(&mult);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&added);
    igraph_vector_int_destroy(&vertex_index);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(proj, &edges, remaining_nodes,
                               /*directed=*/ 0));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_destroy, proj);

    IGRAPH_I_ATTRIBUTE_DESTROY(proj);
    IGRAPH_I_ATTRIBUTE_COPY(proj, graph, 1, 0, 0);
    IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, proj, &vertex_perm));
    igraph_vector_int_destroy(&vertex_perm);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bipartite_projection
 * \brief Create one or both projections of a bipartite (two-mode) network.
 *
 * Creates one or both projections of a bipartite graph.
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
                                igraph_integer_t probe1) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    /* t1 is -1 if proj1 is omitted, it is 0 if it belongs to type zero,
       it is 1 if it belongs to type one. The same for t2 */
    int t1, t2;

    if (igraph_vector_bool_size(types) != no_of_nodes) {
        IGRAPH_ERROR("Invalid bipartite type vector size", IGRAPH_EINVAL);
    }

    if (probe1 >= no_of_nodes) {
        IGRAPH_ERROR("No such vertex to probe", IGRAPH_EINVAL);
    }

    if (probe1 >= 0 && !proj1) {
        IGRAPH_ERROR("`probe1' given, but `proj1' is a null pointer", IGRAPH_EINVAL);
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

    IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj1, t1, multiplicity1));
    IGRAPH_FINALLY(igraph_destroy, proj1);
    IGRAPH_CHECK(igraph_i_bipartite_projection(graph, types, proj2, t2, multiplicity2));

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_full_bipartite
 * \brief Create a full bipartite network.
 *
 * A bipartite network contains two kinds of vertices and connections
 * are only possible between two vertices of different kind. There are
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
 * to them.
 *
 * \param graph Pointer to an igraph_t object, the graph will be
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
 * \sa \ref igraph_full() for non-bipartite full graphs.
 */

igraph_error_t igraph_full_bipartite(igraph_t *graph,
                          igraph_vector_bool_t *types,
                          igraph_integer_t n1, igraph_integer_t n2,
                          igraph_bool_t directed,
                          igraph_neimode_t mode) {

    igraph_integer_t nn1 = n1, nn2 = n2;
    igraph_integer_t no_of_nodes = nn1 + nn2;
    igraph_vector_int_t edges;
    igraph_integer_t no_of_edges;
    igraph_integer_t ptr = 0;
    igraph_integer_t i, j;

    if (!directed) {
        no_of_edges = nn1 * nn2;
    } else if (mode == IGRAPH_OUT || mode == IGRAPH_IN) {
        no_of_edges = nn1 * nn2;
    } else { /* mode==IGRAPH_ALL */
        no_of_edges = nn1 * nn2 * 2;
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);

    if (!directed || mode == IGRAPH_OUT) {

        for (i = 0; i < nn1; i++) {
            for (j = 0; j < nn2; j++) {
                VECTOR(edges)[ptr++] = i;
                VECTOR(edges)[ptr++] = nn1 + j;
            }
        }

    } else if (mode == IGRAPH_IN) {

        for (i = 0; i < nn1; i++) {
            for (j = 0; j < nn2; j++) {
                VECTOR(edges)[ptr++] = nn1 + j;
                VECTOR(edges)[ptr++] = i;
            }
        }

    } else {

        for (i = 0; i < nn1; i++) {
            for (j = 0; j < nn2; j++) {
                VECTOR(edges)[ptr++] = i;
                VECTOR(edges)[ptr++] = nn1 + j;
                VECTOR(edges)[ptr++] = nn1 + j;
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
        for (i = nn1; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = 1;
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
 * \param directed Boolean scalar, whether to create a directed
 *   graph.
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

    igraph_integer_t no_of_nodes = igraph_vector_bool_size(types);
    igraph_integer_t no_of_edges = igraph_vector_int_size(edges);
    igraph_integer_t i;

    if (no_of_edges % 2 != 0) {
        IGRAPH_ERROR("Invalid (odd) edges vector", IGRAPH_EINVEVECTOR);
    }
    no_of_edges /= 2;

    if (! igraph_vector_int_isininterval(edges, 0, no_of_nodes-1)) {
        IGRAPH_ERROR("Invalid (negative or too large) vertex ID", IGRAPH_EINVVID);
    }

    /* Check bipartiteness */
    for (i = 0; i < no_of_edges * 2; i += 2) {
        igraph_integer_t from = VECTOR(*edges)[i];
        igraph_integer_t to = VECTOR(*edges)[i + 1];
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
 * \function igraph_incidence
 * \brief Creates a bipartite graph from an incidence matrix.
 *
 * A bipartite (or two-mode) graph contains two types of vertices and
 * edges always connect vertices of different types. An incidence
 * matrix is an \em n x \em m matrix, \em n and \em m are the number of vertices
 * of the two types, respectively. Nonzero elements in the matrix denote
 * edges between the two corresponding vertices.
 *
 * </para><para>
 * Note that this function can operate in two modes, depending on the
 * \p multiple argument. If it is \c false, then a single edge is
 * created for every non-zero element in the incidence matrix. If \p
 * multiple is \c true, then the matrix elements are rounded up
 * to the closest non-negative integer to get the number of edges to
 * create between a pair of vertices.
 *
 * </para><para>
 * This function does not create multiple edges if \p multiple is
 * \c false, but might create some if it is \c true.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param types Pointer to an initialized boolean vector, or a null
 *   pointer. If not a null pointer, then the vertex types are stored
 *   here. It is resized as needed.
 * \param incidence The incidence matrix.
 * \param directed Gives whether to create an undirected or a directed
 *   graph.
 * \param mode Specifies the direction of the edges in a directed
 *   graph. If \c IGRAPH_OUT, then edges point from vertices
 *   of the first kind (corresponding to rows) to vertices of the
 *   second kind (corresponding to columns); if \c
 *   IGRAPH_IN, then the opposite direction is realized; if \c
 *   IGRAPH_ALL, then mutual edges will be created.
 * \param multiple How to interpret the incidence matrix elements. See
 *   details below.
 * \return Error code.
 *
 * Time complexity: O(n*m), the size of the incidence matrix.
 */

igraph_error_t igraph_incidence(igraph_t *graph, igraph_vector_bool_t *types,
                     const igraph_matrix_t *incidence,
                     igraph_bool_t directed,
                     igraph_neimode_t mode, igraph_bool_t multiple) {

    igraph_integer_t n1 = igraph_matrix_nrow(incidence);
    igraph_integer_t n2 = igraph_matrix_ncol(incidence);
    igraph_integer_t no_of_nodes = n1 + n2;
    igraph_vector_int_t edges;
    igraph_integer_t i, j, k;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    if (n1 > 0 && n2 > 0 && igraph_matrix_min(incidence) < 0) {
        IGRAPH_ERRORF("Incidence matrix elements should be non-negative, found %g.",
                IGRAPH_EINVAL, igraph_matrix_min(incidence));
    }

    if (multiple) {

        for (i = 0; i < n1; i++) {
            for (j = 0; j < n2; j++) {
                igraph_integer_t elem = ceil(MATRIX(*incidence, i, j));
                igraph_integer_t from, to;

                if (!elem) {
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
                    for (k = 0; k < elem; k++) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                    }
                } else {
                    for (k = 0; k < elem; k++) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, from));
                    }
                }
            }
        }

    } else {

        for (i = 0; i < n1; i++) {
            for (j = 0; j < n2; j++) {
                igraph_integer_t from, to;

                if (MATRIX(*incidence, i, j) != 0) {
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
        for (i = n1; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = 1;
        }
    }

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_get_incidence
 * \brief Convert a bipartite graph into an incidence matrix.
 *
 * \param graph The input graph, edge directions are ignored.
 * \param types Boolean vector containing the vertex types. All vertices
 *   in one part of the graph should have type 0, the others type 1.
 * \param res Pointer to an initialized matrix, the result is stored
 *   here. An element of the matrix gives the number of edges
 *   (irrespectively of their direction) between the two corresponding
 *   vertices. The rows will correspond to vertices with type 0,
 *   the columns correspond to vertices with type 1.
 * \param row_ids Pointer to an initialized vector or a null
 *   pointer. If not a null pointer, then the vertex IDs (in the
 *   graph) corresponding to the rows of the result matrix are stored
 *   here.
 * \param col_ids Pointer to an initialized vector or a null
 *   pointer. If not a null pointer, then the vertex IDs corresponding
 *   to the columns of the result matrix are stored here.
 * \return Error code.
 *
 * Time complexity: O(n*m), n and m are number of vertices of the two
 * different kind.
 *
 * \sa \ref igraph_incidence() for the opposite operation.
 */

igraph_error_t igraph_get_incidence(const igraph_t *graph,
                         const igraph_vector_bool_t *types,
                         igraph_matrix_t *res,
                         igraph_vector_int_t *row_ids,
                         igraph_vector_int_t *col_ids) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t n1 = 0, n2 = 0, i;
    igraph_vector_int_t perm;
    igraph_integer_t p1, p2;
    igraph_integer_t ignored_edges = 0;

    if (igraph_vector_bool_size(types) != no_of_nodes) {
        IGRAPH_ERRORF("Vertex type vector size (%" IGRAPH_PRId ") not equal to number of vertices (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_bool_size(types), no_of_nodes);
    }

    for (i = 0; i < no_of_nodes; i++) {
        n1 += VECTOR(*types)[i] == false ? 1 : 0;
    }
    n2 = no_of_nodes - n1;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&perm, no_of_nodes);

    for (i = 0, p1 = 0, p2 = n1; i < no_of_nodes; i++) {
        VECTOR(perm)[i] = VECTOR(*types)[i] ? p2++ : p1++;
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, n1, n2));
    igraph_matrix_null(res);
    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, i);
        igraph_integer_t to = IGRAPH_TO(graph, i);
        igraph_integer_t from2 = VECTOR(perm)[from];
        igraph_integer_t to2 = VECTOR(perm)[to];
        if (VECTOR(*types)[from] == VECTOR(*types)[to]) {
            ignored_edges++;
        } else if (! VECTOR(*types)[from]) {
            MATRIX(*res, from2, to2 - n1) += 1;
        } else {
            MATRIX(*res, to2, from2 - n1) += 1;
        }
    }
    if (ignored_edges) {
            IGRAPH_WARNINGF("%" IGRAPH_PRId " edges running within partitions were ignored.", ignored_edges);
    }

    if (row_ids) {
        IGRAPH_CHECK(igraph_vector_int_resize(row_ids, n1));
    }
    if (col_ids) {
        IGRAPH_CHECK(igraph_vector_int_resize(col_ids, n2));
    }
    if (row_ids || col_ids) {
        for (i = 0; i < no_of_nodes; i++) {
            if (! VECTOR(*types)[i]) {
                if (row_ids) {
                    igraph_integer_t i2 = VECTOR(perm)[i];
                    VECTOR(*row_ids)[i2] = i;
                }
            } else {
                if (col_ids) {
                    igraph_integer_t i2 = VECTOR(perm)[i];
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

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_char_t seen;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t neis;
    igraph_bool_t bi = true;

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&seen, no_of_nodes);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    for (igraph_integer_t i = 0; bi && i < no_of_nodes; i++) {

        if (VECTOR(seen)[i]) {
            continue;
        }

        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, i));
        VECTOR(seen)[i] = 1;

        while (bi && !igraph_dqueue_int_empty(&Q)) {
            igraph_integer_t n, j;
            igraph_integer_t actnode = igraph_dqueue_int_pop(&Q);
            char acttype = VECTOR(seen)[actnode];

            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, IGRAPH_ALL));
            n = igraph_vector_int_size(&neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t nei = VECTOR(neis)[j];
                if (VECTOR(seen)[nei]) {
                    igraph_integer_t neitype = VECTOR(seen)[nei];
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

    if (res) {
        *res = bi;
    }

    if (types && bi) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, no_of_nodes));
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            VECTOR(*types)[i] = VECTOR(seen)[i] - 1;
        }
    }

    igraph_vector_char_destroy(&seen);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_bipartite_game_gnp(igraph_t *graph, igraph_vector_bool_t *types,
                              igraph_integer_t n1, igraph_integer_t n2,
                              igraph_real_t p, igraph_bool_t directed,
                              igraph_neimode_t mode) {

    igraph_vector_int_t edges, s;
    igraph_integer_t i;

    if (p < 0.0 || p > 1.0) {
        IGRAPH_ERROR("Invalid connection probability", IGRAPH_EINVAL);
    }

    if (types) {
        IGRAPH_CHECK(igraph_vector_bool_resize(types, n1 + n2));
        igraph_vector_bool_null(types);
        for (i = n1; i < n1 + n2; i++) {
            VECTOR(*types)[i] = 1;
        }
    }

    if (p == 0 || n1 * n2 < 1) {
        IGRAPH_CHECK(igraph_empty(graph, n1 + n2, directed));
    } else if (p == 1.0) {
        IGRAPH_CHECK(igraph_full_bipartite(graph, types, n1, n2, directed,
                              mode));
    } else {

        igraph_integer_t to, from, slen;
        igraph_real_t n1_real = n1;  /* for divisions below */
        igraph_real_t n2_real = n2;  /* for divisions below */
        igraph_real_t maxedges, last;
        igraph_integer_t maxedges_int;

        if (!directed || mode != IGRAPH_ALL) {
            maxedges = n1_real * n2_real;
        } else {
            maxedges = 2.0 * n1_real * n2_real;
        }

        if (maxedges > IGRAPH_MAX_EXACT_REAL) {
            IGRAPH_ERROR("Too many vertices, overflow in maximum number of edges.", IGRAPH_EOVERFLOW);
        }
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_i_safe_floor(maxedges * p * 1.1, &maxedges_int));
        IGRAPH_CHECK(igraph_vector_int_reserve(&s, maxedges_int));

        RNG_BEGIN();

        last = RNG_GEOM(p);
        while (last < maxedges) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&s, last));
            last += RNG_GEOM(p);
            last += 1;
        }

        RNG_END();

        slen = igraph_vector_int_size(&s);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, slen * 2));

        for (i = 0; i < slen; i++) {
            if (!directed || mode != IGRAPH_ALL) {
                to = floor(VECTOR(s)[i] / n1_real);
                from = VECTOR(s)[i] - to * n1_real;
                to += n1;
            } else {
                igraph_integer_t n1n2 = n1 * n2;
                if (VECTOR(s)[i] < n1n2) {
                    to = floor(VECTOR(s)[i] / n1_real);
                    from = VECTOR(s)[i] - to * n1_real;
                    to += n1;
                } else {
                    to = floor((VECTOR(s)[i] - n1n2) / n2_real);
                    from = VECTOR(s)[i] - n1n2 - to * n2_real;
                    from += n1;
                }
            }

            if (mode != IGRAPH_IN) {
                igraph_vector_int_push_back(&edges, from);
                igraph_vector_int_push_back(&edges, to);
            } else {
                igraph_vector_int_push_back(&edges, to);
                igraph_vector_int_push_back(&edges, from);
            }
        }

        igraph_vector_int_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_CHECK(igraph_create(graph, &edges, n1 + n2, directed));
        igraph_vector_int_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_bipartite_game_gnm(igraph_t *graph, igraph_vector_bool_t *types,
                              igraph_integer_t n1, igraph_integer_t n2,
                              igraph_integer_t m, igraph_bool_t directed,
                              igraph_neimode_t mode) {
    igraph_vector_int_t edges;
    igraph_vector_int_t s;

    if (n1 < 0 || n2 < 0) {
        IGRAPH_ERROR("Invalid number of vertices.", IGRAPH_EINVAL);
    }
    if (m < 0) {
        IGRAPH_ERROR("Invalid number of edges.", IGRAPH_EINVAL);
    }

    if (types) {
        igraph_integer_t i;
        IGRAPH_CHECK(igraph_vector_bool_resize(types, n1 + n2));
        igraph_vector_bool_null(types);
        for (i = n1; i < n1 + n2; i++) {
            VECTOR(*types)[i] = 1;
        }
    }

    if (m == 0 || n1 * n2 == 0) {
        if (m > 0) {
            IGRAPH_ERROR("Too many edges requested compared to the number of vertices.", IGRAPH_EINVAL);
        }
        IGRAPH_CHECK(igraph_empty(graph, n1 + n2, directed));
    } else {
        igraph_integer_t i;
        igraph_real_t maxedges;

        if (!directed || mode != IGRAPH_ALL) {
            maxedges = (igraph_real_t) n1 * (igraph_real_t) n2;
        } else {
            maxedges = 2.0 * (igraph_real_t) n1 * (igraph_real_t) n2;
        }

        if (m > maxedges) {
            IGRAPH_ERROR("Too many edges requested compared to the number of vertices.", IGRAPH_EINVAL);
        }

        if (maxedges == m) {
            IGRAPH_CHECK(igraph_full_bipartite(graph, types, n1, n2,
                                  directed, mode));
        } else {

            igraph_integer_t to, from;
            igraph_real_t n1_real = (igraph_real_t) n1;  /* for divisions below */
            igraph_real_t n2_real = (igraph_real_t) n2;  /* for divisions below */

            IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&s, 0);
            IGRAPH_CHECK(igraph_random_sample(&s, 0, maxedges - 1, m));
            IGRAPH_CHECK(igraph_vector_int_reserve(&edges, igraph_vector_int_size(&s) * 2));

            for (i = 0; i < m; i++) {
                if (!directed || mode != IGRAPH_ALL) {
                    to = floor(VECTOR(s)[i] / n1_real);
                    from = VECTOR(s)[i] - to * n1_real;
                    to += n1;
                } else {
                    igraph_integer_t n1n2 = n1 * n2;
                    if (VECTOR(s)[i] < n1n2) {
                        to = floor(VECTOR(s)[i] / n1_real);
                        from = VECTOR(s)[i] - to * n1_real;
                        to += n1;
                    } else {
                        to = floor((VECTOR(s)[i] - n1n2) / n2_real);
                        from = VECTOR(s)[i] - n1n2 - to * n2_real;
                        from += n1;
                    }
                }

                if (mode != IGRAPH_IN) {
                    igraph_vector_int_push_back(&edges, from);
                    igraph_vector_int_push_back(&edges, to);
                } else {
                    igraph_vector_int_push_back(&edges, to);
                    igraph_vector_int_push_back(&edges, from);
                }
            }

            igraph_vector_int_destroy(&s);
            IGRAPH_FINALLY_CLEAN(1);
            IGRAPH_CHECK(igraph_create(graph, &edges, n1 + n2, directed));
            igraph_vector_int_destroy(&edges);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bipartite_game
 * \brief Generate a bipartite random graph (similar to Erdős-Rényi).
 *
 * Similarly to unipartite (one-mode) networks, we can define the
 * G(n,p), and G(n,m) graph classes for bipartite graphs, via their
 * generating process. In G(n,p) every possible edge between top and
 * bottom vertices is realized with probability p, independently of the
 * rest of the edges. In G(n,m), we uniformly choose m edges to
 * realize.
 *
 * \param graph Pointer to an uninitialized igraph graph, the result
 *    is stored here.
 * \param types Pointer to an initialized boolean vector, or a null
 *    pointer. If not a null pointer, then the vertex types are stored
 *    here. Bottom vertices come first, n1 of them, then n2 top
 *    vertices.
 * \param type The type of the random graph, possible values:
 *        \clist
 *        \cli IGRAPH_ERDOS_RENYI_GNM
 *          G(n,m) graph,
 *          m edges are
 *          selected uniformly randomly in a graph with
 *          n vertices.
 *        \cli IGRAPH_ERDOS_RENYI_GNP
 *          G(n,p) graph,
 *          every possible edge is included in the graph with
 *          probability p.
 *        \endclist
 * \param n1 The number of bottom vertices.
 * \param n2 The number of top vertices.
 * \param p The connection probability for G(n,p) graphs. It is
 *     ignored for G(n,m) graphs.
 * \param m The number of edges for G(n,m) graphs. It is ignored for
 *     G(n,p) graphs.
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
 * \return Error code.
 *
 * \sa \ref igraph_erdos_renyi_game.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

igraph_error_t igraph_bipartite_game(igraph_t *graph, igraph_vector_bool_t *types,
                          igraph_erdos_renyi_t type,
                          igraph_integer_t n1, igraph_integer_t n2,
                          igraph_real_t p, igraph_integer_t m,
                          igraph_bool_t directed, igraph_neimode_t mode) {

    if (n1 < 0 || n2 < 0) {
        IGRAPH_ERROR("Invalid number of vertices for bipartite game.", IGRAPH_EINVAL);
    }

    if (type == IGRAPH_ERDOS_RENYI_GNP) {
        return igraph_bipartite_game_gnp(graph, types, n1, n2, p, directed, mode);
    } else if (type == IGRAPH_ERDOS_RENYI_GNM) {
        return igraph_bipartite_game_gnm(graph, types, n1, n2, m, directed, mode);
    } else {
        IGRAPH_ERROR("Invalid bipartite game type.", IGRAPH_EINVAL);
    }
}
