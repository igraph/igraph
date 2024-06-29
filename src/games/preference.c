/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_constructors.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_vector_list.h"

#include "core/interruption.h"
#include "math/safe_intop.h"

#include <math.h> /* for sqrt and floor */

/**
 * \function igraph_preference_game
 * \brief Generates a graph with vertex types and connection preferences.
 *
 * </para><para>
 * This is practically the nongrowing variant of
 * \ref igraph_establishment_game(). A given number of vertices are
 * generated. Every vertex is assigned to a vertex type according to
 * the given type probabilities. Finally, every
 * vertex pair is evaluated and an edge is created between them with a
 * probability depending on the types of the vertices involved.
 *
 * </para><para>
 * In other words, this function generates a graph according to a
 * block-model. Vertices are divided into groups (or blocks), and
 * the probability the two vertices are connected depends on their
 * groups only.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param type_dist Vector giving the distribution of vertex types. If
 *   \c NULL, all vertex types will have equal probability. See also the
 *   \p fixed_sizes argument.
 * \param fixed_sizes Boolean. If true, then the number of vertices with a
 *   given vertex type is fixed and the \p type_dist argument gives these
 *   numbers for each vertex type. If true, and \p type_dist is \c NULL,
 *   then the function tries to make vertex groups of the same size. If this
 *   is not possible, then some groups will have an extra vertex.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types. This should be symmetric if the requested
 *   graph is undirected.
 * \param node_type_vec A vector where the individual generated vertex types
 *   will be stored. If \c NULL, the vertex types won't be saved.
 * \param directed Logical, whether to generate a directed graph. If undirected
 *   graphs are requested, only the lower left triangle of the preference
 *   matrix is considered.
 * \param loops Logical, whether loop edges are allowed.
 * \return Error code.
 *
 * Added in version 0.3.</para><para>
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_asymmetric_preference_game(),
 * \ref igraph_establishment_game(), \ref igraph_callaway_traits_game()
 */

igraph_error_t igraph_preference_game(igraph_t *graph, igraph_integer_t nodes,
                           igraph_integer_t types,
                           const igraph_vector_t *type_dist,
                           igraph_bool_t fixed_sizes,
                           const igraph_matrix_t *pref_matrix,
                           igraph_vector_int_t *node_type_vec,
                           igraph_bool_t directed,
                           igraph_bool_t loops) {

    igraph_integer_t i, j, no_reserved_edges;
    igraph_vector_int_t edges;
    igraph_vector_t s;
    igraph_vector_int_t* nodetypes;
    igraph_vector_int_list_t vids_by_type;
    igraph_real_t maxcum, maxedges;

    if (nodes < 0) {
        IGRAPH_ERROR("The number of vertices must be non-negative.", IGRAPH_EINVAL);
    }

    if (types < 1) {
        IGRAPH_ERROR("The number of vertex types must be at least 1.", IGRAPH_EINVAL);
    }

    if (type_dist) {
        igraph_real_t lo;

        if (igraph_vector_size(type_dist) != types) {
            IGRAPH_ERROR("The vertex type distribution vector must agree in length with the number of types.",
                         IGRAPH_EINVAL);
        }

        lo = igraph_vector_min(type_dist);
        if (lo < 0) {
            IGRAPH_ERROR("The vertex type distribution vector must not contain negative values.", IGRAPH_EINVAL);
        }
        if (isnan(lo)) {
            IGRAPH_ERROR("The vertex type distribution vector must not contain NaN.", IGRAPH_EINVAL);
        }
    }

    if (igraph_matrix_nrow(pref_matrix) != types || igraph_matrix_ncol(pref_matrix) != types) {
        IGRAPH_ERROR("The preference matrix must be square and agree in dimensions with the number of types.", IGRAPH_EINVAL);
    }

    {
        igraph_real_t lo, hi;
        igraph_matrix_minmax(pref_matrix, &lo, &hi); /* matrix size is at least 1x1, safe to call minmax */

        if (lo < 0 || hi > 1) {
            IGRAPH_ERROR("The preference matrix must contain probabilities in [0, 1].", IGRAPH_EINVAL);
        }
        if (isnan(lo) || isnan(hi)) {
            IGRAPH_ERROR("The preference matrix must not contain NaN.", IGRAPH_EINVAL);
        }
    }

    if (! directed && ! igraph_matrix_is_symmetric(pref_matrix)) {
        IGRAPH_ERROR("The preference matrix must be symmetric when generating undirected graphs.", IGRAPH_EINVAL);
    }

    if (fixed_sizes && type_dist) {
        if (igraph_vector_sum(type_dist) != nodes) {
            IGRAPH_ERROR("Invalid group sizes, their sum must match the number of vertices.", IGRAPH_EINVAL);
        }
    }

    if (node_type_vec) {
        IGRAPH_CHECK(igraph_vector_int_resize(node_type_vec, nodes));
        nodetypes = node_type_vec;
    } else {
        nodetypes = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(nodetypes, "Insufficient memory for preference_game.");
        IGRAPH_FINALLY(igraph_free, nodetypes);
        IGRAPH_VECTOR_INT_INIT_FINALLY(nodetypes, nodes);
    }

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&vids_by_type, types);

    RNG_BEGIN();

    if (!fixed_sizes) {

        igraph_vector_t cumdist;
        IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types + 1);

        VECTOR(cumdist)[0] = 0;
        if (type_dist) {
            for (i = 0; i < types; i++) {
                VECTOR(cumdist)[i + 1] = VECTOR(cumdist)[i] + VECTOR(*type_dist)[i];
            }
        } else {
            for (i = 0; i < types; i++) {
                VECTOR(cumdist)[i + 1] = i + 1;
            }
        }
        maxcum = igraph_vector_tail(&cumdist);

        for (i = 0; i < nodes; i++) {
            igraph_integer_t type1;
            igraph_real_t uni1 = RNG_UNIF(0, maxcum);
            igraph_vector_binsearch(&cumdist, uni1, &type1);
            VECTOR(*nodetypes)[i] = type1 - 1;
            IGRAPH_CHECK(igraph_vector_int_push_back(
                igraph_vector_int_list_get_ptr(&vids_by_type, type1 - 1), i
            ));
        }

        igraph_vector_destroy(&cumdist);
        IGRAPH_FINALLY_CLEAN(1);

    } else {
        igraph_integer_t an = 0;
        if (type_dist) {
            for (i = 0; i < types; i++) {
                igraph_integer_t no = VECTOR(*type_dist)[i];
                igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&vids_by_type, i);
                for (j = 0; j < no && an < nodes; j++) {
                    VECTOR(*nodetypes)[an] = i;
                    IGRAPH_CHECK(igraph_vector_int_push_back(v, an));
                    an++;
                }
            }
        } else {
            igraph_integer_t size_of_one_group = nodes / types;
            igraph_integer_t num_groups_with_one_extra_node = nodes - size_of_one_group * types;
            for (i = 0; i < types; i++) {
                igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&vids_by_type, i);
                for (j = 0; j < size_of_one_group; j++) {
                    VECTOR(*nodetypes)[an] = i;
                    IGRAPH_CHECK(igraph_vector_int_push_back(v, an));
                    an++;
                }
                if (i < num_groups_with_one_extra_node) {
                    VECTOR(*nodetypes)[an] = i;
                    IGRAPH_CHECK(igraph_vector_int_push_back(v, an));
                    an++;
                }
            }
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);

    for (i = 0; i < types; i++) {
        for (j = 0; j < types; j++) {
            /* Generating the random subgraph between vertices of type i and j */
            igraph_integer_t k, l, l_x2;
            igraph_real_t p, last;
            igraph_vector_int_t *v1, *v2;
            igraph_integer_t v1_size, v2_size;

            IGRAPH_ALLOW_INTERRUPTION();

            v1 = igraph_vector_int_list_get_ptr(&vids_by_type, i);
            v2 = igraph_vector_int_list_get_ptr(&vids_by_type, j);
            v1_size = igraph_vector_int_size(v1);
            v2_size = igraph_vector_int_size(v2);

            p = MATRIX(*pref_matrix, i, j);
            igraph_vector_clear(&s);
            if (i != j) {
                /* The two vertex sets are disjoint, this is the easier case */
                if (i > j && !directed) {
                    continue;
                }
                maxedges = ((igraph_real_t) v1_size) * v2_size;
            } else {
                if (directed && loops) {
                    maxedges = ((igraph_real_t) v1_size) * v1_size;
                } else if (directed && !loops) {
                    maxedges = ((igraph_real_t) v1_size) * (v1_size - 1);
                } else if (!directed && loops) {
                    maxedges = ((igraph_real_t) v1_size) * (v1_size + 1) / 2;
                } else {
                    maxedges = ((igraph_real_t) v1_size) * (v1_size - 1) / 2;
                }
            }

            if (maxedges > IGRAPH_MAX_EXACT_REAL) {
                IGRAPH_ERROR("Too many vertices, overflow in maximum number of edges.", IGRAPH_EOVERFLOW);
            }

            IGRAPH_CHECK(igraph_i_safe_floor(maxedges * p * 1.1, &no_reserved_edges));
            IGRAPH_CHECK(igraph_vector_reserve(&s, no_reserved_edges));

            last = RNG_GEOM(p);
            while (last < maxedges) {
                IGRAPH_CHECK(igraph_vector_push_back(&s, last));
                last += RNG_GEOM(p);
                last += 1;
            }
            l = igraph_vector_size(&s);

            IGRAPH_SAFE_MULT(l, 2, &l_x2);
            IGRAPH_SAFE_ADD(igraph_vector_int_size(&edges), l_x2, &no_reserved_edges);
            IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_reserved_edges));

            if (i != j) {
                /* Generating the subgraph between vertices of type i and j */
                for (k = 0; k < l; k++) {
                    igraph_integer_t to = floor(VECTOR(s)[k] / v1_size);
                    igraph_integer_t from = (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                    igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_int_push_back(&edges, VECTOR(*v2)[to]);
                }
            } else {
                /* Generating the subgraph among vertices of type i */
                if (directed && loops) {
                    for (k = 0; k < l; k++) {
                        igraph_integer_t to = floor(VECTOR(s)[k] / v1_size);
                        igraph_integer_t from = (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else if (directed && !loops) {
                    for (k = 0; k < l; k++) {
                        igraph_integer_t to = floor(VECTOR(s)[k] / v1_size);
                        igraph_integer_t from = (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                        if (from == to) {
                            to = v1_size - 1;
                        }
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else if (!directed && loops) {
                    for (k = 0; k < l; k++) {
                        igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[k] + 1) - 1) / 2);
                        igraph_integer_t from = (VECTOR(s)[k] - (((igraph_real_t)to) * (to + 1)) / 2);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else {
                    for (k = 0; k < l; k++) {
                        igraph_integer_t to = floor((sqrt(8 * VECTOR(s)[k] + 1) + 1) / 2);
                        igraph_integer_t from = (VECTOR(s)[k] - (((igraph_real_t)to) * (to - 1)) / 2);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_int_push_back(&edges, VECTOR(*v1)[to]);
                    }
                }
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&s);
    igraph_vector_int_list_destroy(&vids_by_type);
    IGRAPH_FINALLY_CLEAN(2);

    if (node_type_vec == 0) {
        igraph_vector_int_destroy(nodetypes);
        IGRAPH_FREE(nodetypes);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_asymmetric_preference_game
 * \brief Generates a graph with asymmetric vertex types and connection preferences.
 *
 * </para><para>
 * This is the asymmetric variant of \ref igraph_preference_game().
 * A given number of vertices are generated. Every vertex is assigned to an
 * "outgoing" and an "incoming " vertex type according to the given joint
 * type probabilities. Finally, every vertex pair is evaluated and a
 * directed edge is created between them with a probability depending on the
 * "outgoing" type of the source vertex and the "incoming" type of the target
 * vertex.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param no_out_types The number of vertex out-types.
 * \param no_in_types The number of vertex in-types.
 * \param type_dist_matrix Matrix of size <code>out_types * in_types</code>,
 *   giving the joint distribution of vertex types.
 *   If \c NULL, incoming and outgoing vertex types are independent and uniformly
 *   distributed.
 * \param pref_matrix Matrix of size <code>out_types * in_types</code>,
 *   giving the connection probabilities for different vertex types.
 * \param node_type_out_vec A vector where the individual generated "outgoing"
 *   vertex types will be stored. If \c NULL, the vertex types won't be saved.
 * \param node_type_in_vec A vector where the individual generated "incoming"
 *   vertex types will be stored. If \c NULL, the vertex types won't be saved.
 * \param loops Logical, whether loop edges are allowed.
 * \return Error code.
 *
 * Added in version 0.3.</para><para>
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_preference_game()
 */

igraph_error_t igraph_asymmetric_preference_game(igraph_t *graph, igraph_integer_t nodes,
                                      igraph_integer_t no_out_types,
                                      igraph_integer_t no_in_types,
                                      const igraph_matrix_t *type_dist_matrix,
                                      const igraph_matrix_t *pref_matrix,
                                      igraph_vector_int_t *node_type_out_vec,
                                      igraph_vector_int_t *node_type_in_vec,
                                      igraph_bool_t loops) {

    igraph_integer_t i, j, k, no_reserved_edges;
    igraph_vector_int_t edges;
    igraph_vector_t s;
    igraph_vector_t cumdist;
    igraph_vector_int_t intersect;
    igraph_vector_int_t *nodetypes_in;
    igraph_vector_int_t *nodetypes_out;
    igraph_vector_int_list_t vids_by_intype, vids_by_outtype;
    igraph_real_t maxcum, maxedges;

    if (nodes < 0) {
        IGRAPH_ERROR("The number of vertices must not be negative.", IGRAPH_EINVAL);
    }

    if (no_in_types < 1) {
        IGRAPH_ERROR("The number of vertex in-types must be at least 1.", IGRAPH_EINVAL);
    }

    if (no_out_types < 1) {
        IGRAPH_ERROR("The number of vertex out-types must be at least 1.", IGRAPH_EINVAL);
    }

    if (type_dist_matrix) {
        igraph_real_t lo;

        if (igraph_matrix_nrow(type_dist_matrix) != no_out_types ||
            igraph_matrix_ncol(type_dist_matrix) != no_in_types) {
            IGRAPH_ERROR("The type distribution matrix must have dimensions out_types * in_types.", IGRAPH_EINVAL);
        }

        lo = igraph_matrix_min(type_dist_matrix);
        if (lo < 0) {
            IGRAPH_ERROR("The type distribution matrix must not contain negative values.", IGRAPH_EINVAL);
        }
        if (isnan(lo)) {
            IGRAPH_ERROR("The type distribution matrix must not contain NaN.", IGRAPH_EINVAL);
        }
    }

    if (igraph_matrix_nrow(pref_matrix) != no_out_types ||
        igraph_matrix_ncol(pref_matrix) != no_in_types) {
        IGRAPH_ERROR("The preference matrix must have dimensions out_types * in_types.", IGRAPH_EINVAL);
    }

    {
        igraph_real_t lo, hi;
        igraph_matrix_minmax(pref_matrix, &lo, &hi); /* matrix size is at least 1x1, safe to call minmax */

        if (lo < 0 || hi > 1) {
            IGRAPH_ERROR("The preference matrix must contain probabilities in [0, 1].", IGRAPH_EINVAL);
        }
        if (isnan(lo) || isnan(hi)) {
            IGRAPH_ERROR("The preference matrix must not contain NaN.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, no_in_types * no_out_types + 1);

    if (node_type_in_vec) {
        nodetypes_in = node_type_in_vec;
        IGRAPH_CHECK(igraph_vector_int_resize(nodetypes_in, nodes));
    } else {
        nodetypes_in = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(nodetypes_in, "Insufficient memory for asymmetric preference game.");
        IGRAPH_FINALLY(igraph_free, &nodetypes_in);
        IGRAPH_VECTOR_INT_INIT_FINALLY(nodetypes_in, nodes);
    }

    if (node_type_out_vec) {
        nodetypes_out = node_type_out_vec;
        IGRAPH_CHECK(igraph_vector_int_resize(nodetypes_out, nodes));
    } else {
        nodetypes_out = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(nodetypes_out, "Insufficient memory for asymmetric preference game.");
        IGRAPH_FINALLY(igraph_free, &nodetypes_out);
        IGRAPH_VECTOR_INT_INIT_FINALLY(nodetypes_out, nodes);
    }

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&vids_by_intype, no_in_types);
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&vids_by_outtype, no_out_types);

    VECTOR(cumdist)[0] = 0;
    k = 0;
    if (type_dist_matrix) {
        for (j = 0; j < no_in_types; j++) {
            for (i = 0; i < no_out_types; i++) {
                VECTOR(cumdist)[k + 1] = VECTOR(cumdist)[k] + MATRIX(*type_dist_matrix, i, j);
                k++;
            }
        }
    } else {
        for (i = 0; i < no_out_types * no_in_types; i++) {
            VECTOR(cumdist)[i + 1] = i + 1;
        }
    }
    maxcum = igraph_vector_tail(&cumdist);

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        igraph_integer_t in_type, out_type;
        igraph_real_t uni1 = RNG_UNIF(0, maxcum);
        igraph_vector_binsearch(&cumdist, uni1, &in_type);
        out_type = (in_type - 1) % no_out_types;
        in_type = (in_type - 1) / no_out_types;
        VECTOR(*nodetypes_in)[i] = in_type;
        VECTOR(*nodetypes_out)[i] = out_type;
        IGRAPH_CHECK(igraph_vector_int_push_back(
            igraph_vector_int_list_get_ptr(&vids_by_intype, in_type), i
        ));
        IGRAPH_CHECK(igraph_vector_int_push_back(
            igraph_vector_int_list_get_ptr(&vids_by_outtype, out_type), i
        ));
    }

    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&intersect, 0);
    for (i = 0; i < no_out_types; i++) {
        for (j = 0; j < no_in_types; j++) {
            igraph_integer_t kk, l, l_x2;
            igraph_integer_t c = 0;
            igraph_real_t p, last;
            igraph_vector_int_t *v1, *v2;
            igraph_integer_t v1_size, v2_size;

            IGRAPH_ALLOW_INTERRUPTION();

            v1 = igraph_vector_int_list_get_ptr(&vids_by_outtype, i);
            v2 = igraph_vector_int_list_get_ptr(&vids_by_intype, j);
            v1_size = igraph_vector_int_size(v1);
            v2_size = igraph_vector_int_size(v2);

            maxedges = ((igraph_real_t) v1_size) * v2_size;

            if (maxedges > IGRAPH_MAX_EXACT_REAL) {
                IGRAPH_ERROR("Too many vertices, overflow in maximum number of edges.", IGRAPH_EOVERFLOW);
            }

            if (!loops) {
                IGRAPH_CHECK(igraph_vector_int_intersect_sorted(v1, v2, &intersect));
                c = igraph_vector_int_size(&intersect);
                maxedges -= c;
            }

            p = MATRIX(*pref_matrix, i, j);
            igraph_vector_clear(&s);

            IGRAPH_CHECK(igraph_i_safe_floor(maxedges * p * 1.1, &no_reserved_edges));
            IGRAPH_CHECK(igraph_vector_reserve(&s, no_reserved_edges));

            last = RNG_GEOM(p);
            while (last < maxedges) {
                IGRAPH_CHECK(igraph_vector_push_back(&s, last));
                last += RNG_GEOM(p);
                last += 1;
            }
            l = igraph_vector_size(&s);

            IGRAPH_SAFE_MULT(l, 2, &l_x2);
            IGRAPH_SAFE_ADD(igraph_vector_int_size(&edges), l_x2, &no_reserved_edges);
            IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_reserved_edges));


            if (!loops && c > 0) {
                for (kk = 0; kk < l; kk++) {
                    igraph_integer_t to = floor(VECTOR(s)[kk] / v1_size);
                    igraph_integer_t from = (VECTOR(s)[kk] - ((igraph_real_t) to) * v1_size);
                    if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                        /* remap loop edges */
                        to = v2_size - 1;
                        igraph_vector_int_binsearch(&intersect, VECTOR(*v1)[from], &c);
                        from = v1_size - 1;
                        if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                            from--;
                        }
                        while (c > 0) {
                            c--; from--;
                            if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                                from--;
                            }
                        }
                    }
                    igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_int_push_back(&edges, VECTOR(*v2)[to]);
                }
            } else {
                for (kk = 0; kk < l; kk++) {
                    igraph_integer_t to = floor(VECTOR(s)[kk] / v1_size);
                    igraph_integer_t from = (VECTOR(s)[kk] - ((igraph_real_t)to) * v1_size);
                    igraph_vector_int_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_int_push_back(&edges, VECTOR(*v2)[to]);
                }
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&s);
    igraph_vector_int_destroy(&intersect);
    igraph_vector_int_list_destroy(&vids_by_intype);
    igraph_vector_int_list_destroy(&vids_by_outtype);
    IGRAPH_FINALLY_CLEAN(4);

    if (node_type_out_vec == 0) {
        igraph_vector_int_destroy(nodetypes_out);
        IGRAPH_FREE(nodetypes_out);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (node_type_in_vec == 0) {
        igraph_vector_int_destroy(nodetypes_in);
        IGRAPH_FREE(nodetypes_in);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, 1));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
