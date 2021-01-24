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

#include "core/interruption.h"

static void igraph_i_preference_game_free_vids_by_type(igraph_vector_ptr_t *vecs) {
    igraph_long_t i = 0, n;
    igraph_vector_t *v;

    n = (igraph_long_t) igraph_vector_ptr_size(vecs);
    for (i = 0; i < n; i++) {
        v = (igraph_vector_t*)VECTOR(*vecs)[i];
        if (v) {
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_destroy_all(vecs);
}

/**
 * \function igraph_preference_game
 * \brief Generates a graph with vertex types and connection preferences
 *
 * </para><para>
 * This is practically the nongrowing variant of \ref
 * igraph_establishment_game. A given number of vertices are
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
 *   \c fixed_sizes argument.
 * \param fixed_sizes Boolean. If true, then the number of vertices with a
 *   given vertex type is fixed and the \c type_dist argument gives these
 *   numbers for each vertex type. If true, and \c type_dist is \c NULL,
 *   then the function tries to make vertex groups of the same size. If this
 *   is not possible, then some groups will have an extra vertex.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types. This should be symmetric if the requested
 *   graph is undirected.
 * \param node_type_vec A vector where the individual generated vertex types
 *   will be stored. If \c NULL , the vertex types won't be saved.
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
 * \sa igraph_establishment_game()
 */

igraph_long_t igraph_preference_game(igraph_t *graph, igraph_long_t nodes,
                           igraph_long_t types,
                           const igraph_vector_t *type_dist,
                           igraph_bool_t fixed_sizes,
                           const igraph_matrix_t *pref_matrix,
                           igraph_vector_t *node_type_vec,
                           igraph_bool_t directed,
                           igraph_bool_t loops) {

    igraph_long_t i, j;
    igraph_vector_t edges, s;
    igraph_vector_t* nodetypes;
    igraph_vector_ptr_t vids_by_type;
    igraph_real_t maxcum, maxedges;

    if (types < 1) {
        IGRAPH_ERROR("types must be >= 1", IGRAPH_EINVAL);
    }
    if (nodes < 0) {
        IGRAPH_ERROR("nodes must be >= 0", IGRAPH_EINVAL);
    }
    if (type_dist && igraph_vector_size(type_dist) != types) {
        if (igraph_vector_size(type_dist) > types) {
            IGRAPH_WARNING("length of type_dist > types, type_dist will be trimmed");
        } else {
            IGRAPH_ERROR("type_dist vector too short", IGRAPH_EINVAL);
        }
    }
    if (igraph_matrix_nrow(pref_matrix) < types ||
        igraph_matrix_ncol(pref_matrix) < types) {
        IGRAPH_ERROR("pref_matrix too small", IGRAPH_EINVAL);
    }

    if (fixed_sizes && type_dist) {
        if (igraph_vector_sum(type_dist) != nodes) {
            IGRAPH_ERROR("Invalid group sizes, their sum must match the number"
                         " of vertices", IGRAPH_EINVAL);
        }
    }

    if (node_type_vec) {
        IGRAPH_CHECK(igraph_vector_resize(node_type_vec, nodes));
        nodetypes = node_type_vec;
    } else {
        nodetypes = igraph_Calloc(1, igraph_vector_t);
        if (nodetypes == 0) {
            IGRAPH_ERROR("preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, nodetypes);
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes, nodes);
    }

    IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_type, types));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_type);
    for (i = 0; i < types; i++) {
        VECTOR(vids_by_type)[i] = igraph_Calloc(1, igraph_vector_t);
        if (VECTOR(vids_by_type)[i] == 0) {
            IGRAPH_ERROR("preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_type)[i], 0));
    }
    IGRAPH_FINALLY_CLEAN(1);   /* removing igraph_vector_ptr_destroy_all */
    IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_type);

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
            igraph_long_t type1;
            igraph_real_t uni1 = RNG_UNIF(0, maxcum);
            igraph_vector_binsearch(&cumdist, uni1, &type1);
            VECTOR(*nodetypes)[i] = type1 - 1;
            IGRAPH_CHECK(igraph_vector_push_back(
                             (igraph_vector_t*)VECTOR(vids_by_type)[type1 - 1], i));
        }

        igraph_vector_destroy(&cumdist);
        IGRAPH_FINALLY_CLEAN(1);

    } else {

        igraph_long_t an = 0;
        if (type_dist) {
            for (i = 0; i < types; i++) {
                igraph_long_t no = (igraph_long_t) VECTOR(*type_dist)[i];
                igraph_vector_t *v = VECTOR(vids_by_type)[i];
                for (j = 0; j < no && an < nodes; j++) {
                    VECTOR(*nodetypes)[an] = i;
                    IGRAPH_CHECK(igraph_vector_push_back(v, an));
                    an++;
                }
            }
        } else {
            igraph_long_t fixno = (igraph_long_t) ceil( (double)nodes / types);
            for (i = 0; i < types; i++) {
                igraph_vector_t *v = VECTOR(vids_by_type)[i];
                for (j = 0; j < fixno && an < nodes; j++) {
                    VECTOR(*nodetypes)[an++] = i;
                    IGRAPH_CHECK(igraph_vector_push_back(v, an));
                    an++;
                }
            }
        }

    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);

    for (i = 0; i < types; i++) {
        for (j = 0; j < types; j++) {
            /* Generating the random subgraph between vertices of type i and j */
            igraph_long_t k, l;
            igraph_real_t p, last;
            igraph_vector_t *v1, *v2;
            igraph_long_t v1_size, v2_size;

            IGRAPH_ALLOW_INTERRUPTION();

            v1 = (igraph_vector_t*)VECTOR(vids_by_type)[i];
            v2 = (igraph_vector_t*)VECTOR(vids_by_type)[j];
            v1_size = igraph_vector_size(v1);
            v2_size = igraph_vector_size(v2);

            p = MATRIX(*pref_matrix, i, j);
            igraph_vector_clear(&s);
            if (i != j) {
                /* The two vertex sets are disjoint, this is the easier case */
                if (i > j && !directed) {
                    continue;
                }
                maxedges = v1_size * v2_size;
            } else {
                if (directed && loops) {
                    maxedges = v1_size * v1_size;
                } else if (directed && !loops) {
                    maxedges = v1_size * (v1_size - 1);
                } else if (!directed && loops) {
                    maxedges = v1_size * (v1_size + 1) / 2;
                } else {
                    maxedges = v1_size * (v1_size - 1) / 2;
                }
            }

            IGRAPH_CHECK(igraph_vector_reserve(&s, (igraph_long_t) (maxedges * p * 1.1)));

            last = RNG_GEOM(p);
            while (last < maxedges) {
                IGRAPH_CHECK(igraph_vector_push_back(&s, last));
                last += RNG_GEOM(p);
                last += 1;
            }
            l = igraph_vector_size(&s);

            IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&edges) + l * 2));

            if (i != j) {
                /* Generating the subgraph between vertices of type i and j */
                for (k = 0; k < l; k++) {
                    igraph_long_t to = (igraph_long_t) floor(VECTOR(s)[k] / v1_size);
                    igraph_long_t from = (igraph_long_t) (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                    igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
                }
            } else {
                /* Generating the subgraph among vertices of type i */
                if (directed && loops) {
                    for (k = 0; k < l; k++) {
                        igraph_long_t to = (igraph_long_t) floor(VECTOR(s)[k] / v1_size);
                        igraph_long_t from = (igraph_long_t) (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else if (directed && !loops) {
                    for (k = 0; k < l; k++) {
                        igraph_long_t to = (igraph_long_t) floor(VECTOR(s)[k] / v1_size);
                        igraph_long_t from = (igraph_long_t) (VECTOR(s)[k] - ((igraph_real_t)to) * v1_size);
                        if (from == to) {
                            to = v1_size - 1;
                        }
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else if (!directed && loops) {
                    for (k = 0; k < l; k++) {
                        igraph_long_t to = (igraph_long_t) floor((sqrt(8 * VECTOR(s)[k] + 1) - 1) / 2);
                        igraph_long_t from = (igraph_long_t) (VECTOR(s)[k] - (((igraph_real_t)to) * (to + 1)) / 2);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                } else {
                    for (k = 0; k < l; k++) {
                        igraph_long_t to = (igraph_long_t) floor((sqrt(8 * VECTOR(s)[k] + 1) + 1) / 2);
                        igraph_long_t from = (igraph_long_t) (VECTOR(s)[k] - (((igraph_real_t)to) * (to - 1)) / 2);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                        igraph_vector_push_back(&edges, VECTOR(*v1)[to]);
                    }
                }
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&s);
    igraph_i_preference_game_free_vids_by_type(&vids_by_type);
    IGRAPH_FINALLY_CLEAN(2);

    if (node_type_vec == 0) {
        igraph_vector_destroy(nodetypes);
        igraph_Free(nodetypes);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_asymmetric_preference_game
 * \brief Generates a graph with asymmetric vertex types and connection preferences
 *
 * </para><para>
 * This is the asymmetric variant of \ref igraph_preference_game() .
 * A given number of vertices are generated. Every vertex is assigned to an
 * "incoming" and an "outgoing" vertex type according to the given joint
 * type probabilities. Finally, every vertex pair is evaluated and a
 * directed edge is created between them with a probability depending on the
 * "outgoing" type of the source vertex and the "incoming" type of the target
 * vertex.
 *
 * \param graph Pointer to an uninitialized graph.
 * \param nodes The number of vertices in the graph.
 * \param types The number of vertex types.
 * \param type_dist_matrix Matrix giving the joint distribution of vertex types.
 *   If null, incoming and outgoing vertex types are independent and uniformly
 *   distributed.
 * \param pref_matrix Matrix giving the connection probabilities for
 *   different vertex types.
 * \param node_type_in_vec A vector where the individual generated "incoming"
 *   vertex types will be stored. If NULL, the vertex types won't be saved.
 * \param node_type_out_vec A vector where the individual generated "outgoing"
 *   vertex types will be stored. If NULL, the vertex types won't be saved.
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

igraph_long_t igraph_asymmetric_preference_game(igraph_t *graph, igraph_long_t nodes,
                                      igraph_long_t types,
                                      igraph_matrix_t *type_dist_matrix,
                                      igraph_matrix_t *pref_matrix,
                                      igraph_vector_t *node_type_in_vec,
                                      igraph_vector_t *node_type_out_vec,
                                      igraph_bool_t loops) {

    igraph_long_t i, j, k;
    igraph_vector_t edges, cumdist, s, intersect;
    igraph_vector_t *nodetypes_in;
    igraph_vector_t *nodetypes_out;
    igraph_vector_ptr_t vids_by_intype, vids_by_outtype;
    igraph_real_t maxcum, maxedges;

    if (types < 1) {
        IGRAPH_ERROR("types must be >= 1", IGRAPH_EINVAL);
    }
    if (nodes < 0) {
        IGRAPH_ERROR("nodes must be >= 0", IGRAPH_EINVAL);
    }
    if (type_dist_matrix) {
        if (igraph_matrix_nrow(type_dist_matrix) < types ||
            igraph_matrix_ncol(type_dist_matrix) < types) {
            IGRAPH_ERROR("type_dist_matrix too small", IGRAPH_EINVAL);
        } else if (igraph_matrix_nrow(type_dist_matrix) > types ||
                   igraph_matrix_ncol(type_dist_matrix) > types) {
            IGRAPH_WARNING("type_dist_matrix will be trimmed");
        }
    }
    if (igraph_matrix_nrow(pref_matrix) < types ||
        igraph_matrix_ncol(pref_matrix) < types) {
        IGRAPH_ERROR("pref_matrix too small", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&cumdist, types * types + 1);

    if (node_type_in_vec) {
        nodetypes_in = node_type_in_vec;
        IGRAPH_CHECK(igraph_vector_resize(nodetypes_in, nodes));
    } else {
        nodetypes_in = igraph_Calloc(1, igraph_vector_t);
        if (nodetypes_in == 0) {
            IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes_in, nodes);
    }

    if (node_type_out_vec) {
        nodetypes_out = node_type_out_vec;
        IGRAPH_CHECK(igraph_vector_resize(nodetypes_out, nodes));
    } else {
        nodetypes_out = igraph_Calloc(1, igraph_vector_t);
        if (nodetypes_out == 0) {
            IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_VECTOR_INIT_FINALLY(nodetypes_out, nodes);
    }

    IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_intype, types));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_intype);
    IGRAPH_CHECK(igraph_vector_ptr_init(&vids_by_outtype, types));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &vids_by_outtype);
    for (i = 0; i < types; i++) {
        VECTOR(vids_by_intype)[i] = igraph_Calloc(1, igraph_vector_t);
        VECTOR(vids_by_outtype)[i] = igraph_Calloc(1, igraph_vector_t);
        if (VECTOR(vids_by_intype)[i] == 0 || VECTOR(vids_by_outtype)[i] == 0) {
            IGRAPH_ERROR("asymmetric_preference_game failed", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_intype)[i], 0));
        IGRAPH_CHECK(igraph_vector_init(VECTOR(vids_by_outtype)[i], 0));
    }
    IGRAPH_FINALLY_CLEAN(2);   /* removing igraph_vector_ptr_destroy_all */
    IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_intype);
    IGRAPH_FINALLY(igraph_i_preference_game_free_vids_by_type, &vids_by_outtype);

    VECTOR(cumdist)[0] = 0;
    if (type_dist_matrix) {
        for (i = 0, k = 0; i < types; i++) {
            for (j = 0; j < types; j++, k++) {
                VECTOR(cumdist)[k + 1] = VECTOR(cumdist)[k] + MATRIX(*type_dist_matrix, i, j);
            }
        }
    } else {
        for (i = 0; i < types * types; i++) {
            VECTOR(cumdist)[i + 1] = i + 1;
        }
    }
    maxcum = igraph_vector_tail(&cumdist);

    RNG_BEGIN();

    for (i = 0; i < nodes; i++) {
        igraph_long_t type1, type2;
        igraph_real_t uni1 = RNG_UNIF(0, maxcum);
        igraph_vector_binsearch(&cumdist, uni1, &type1);
        type2 = (type1 - 1) % (igraph_long_t)types;
        type1 = (type1 - 1) / (igraph_long_t)types;
        VECTOR(*nodetypes_in)[i] = type1;
        VECTOR(*nodetypes_out)[i] = type2;
        IGRAPH_CHECK(igraph_vector_push_back(
                         (igraph_vector_t*)VECTOR(vids_by_intype)[type1], i));
        IGRAPH_CHECK(igraph_vector_push_back(
                         (igraph_vector_t*)VECTOR(vids_by_outtype)[type2], i));
    }

    igraph_vector_destroy(&cumdist);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&intersect, 0);
    for (i = 0; i < types; i++) {
        for (j = 0; j < types; j++) {
            igraph_long_t kk, l, c = 0;
            igraph_real_t p, last;
            igraph_vector_t *v1, *v2;
            igraph_long_t v1_size, v2_size;

            IGRAPH_ALLOW_INTERRUPTION();

            v1 = (igraph_vector_t*)VECTOR(vids_by_outtype)[i];
            v2 = (igraph_vector_t*)VECTOR(vids_by_intype)[j];
            v1_size = igraph_vector_size(v1);
            v2_size = igraph_vector_size(v2);

            maxedges = v1_size * v2_size;
            if (!loops) {
                IGRAPH_CHECK(igraph_vector_intersect_sorted(v1, v2, &intersect));
                c = igraph_vector_size(&intersect);
                maxedges -= c;
            }

            p = MATRIX(*pref_matrix, i, j);
            igraph_vector_clear(&s);
            IGRAPH_CHECK(igraph_vector_reserve(&s, (igraph_long_t) (maxedges * p * 1.1)));

            last = RNG_GEOM(p);
            while (last < maxedges) {
                IGRAPH_CHECK(igraph_vector_push_back(&s, last));
                last += RNG_GEOM(p);
                last += 1;
            }
            l = igraph_vector_size(&s);

            IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&edges) + l * 2));

            if (!loops && c > 0) {
                for (kk = 0; kk < l; kk++) {
                    igraph_long_t to = (igraph_long_t) floor(VECTOR(s)[kk] / v1_size);
                    igraph_long_t from = (igraph_long_t) (VECTOR(s)[kk] - ((igraph_real_t)to) * v1_size);
                    if (VECTOR(*v1)[from] == VECTOR(*v2)[to]) {
                        /* remap loop edges */
                        to = v2_size - 1;
                        igraph_vector_binsearch(&intersect, VECTOR(*v1)[from], &c);
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
                    igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
                }
            } else {
                for (kk = 0; kk < l; kk++) {
                    igraph_long_t to = (igraph_long_t) floor(VECTOR(s)[kk] / v1_size);
                    igraph_long_t from = (igraph_long_t) (VECTOR(s)[kk] - ((igraph_real_t)to) * v1_size);
                    igraph_vector_push_back(&edges, VECTOR(*v1)[from]);
                    igraph_vector_push_back(&edges, VECTOR(*v2)[to]);
                }
            }
        }
    }

    RNG_END();

    igraph_vector_destroy(&s);
    igraph_vector_destroy(&intersect);
    igraph_i_preference_game_free_vids_by_type(&vids_by_intype);
    igraph_i_preference_game_free_vids_by_type(&vids_by_outtype);
    IGRAPH_FINALLY_CLEAN(4);

    if (node_type_out_vec == 0) {
        igraph_vector_destroy(nodetypes_out);
        igraph_Free(nodetypes_out);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (node_type_in_vec == 0) {
        igraph_vector_destroy(nodetypes_in);
        igraph_Free(nodetypes_in);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, 1));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
