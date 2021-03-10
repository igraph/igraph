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
#include "igraph_psumtree.h"
#include "igraph_random.h"
#include "igraph_interface.h"

typedef struct {
    long int no;
    igraph_psumtree_t *sumtrees;
} igraph_i_citing_cited_type_game_struct_t;

static void igraph_i_citing_cited_type_game_free (
        igraph_i_citing_cited_type_game_struct_t *s);

/**
 * \function igraph_lastcit_game
 * \brief Simulates a citation network, based on time passed since the last citation.
 *
 * This is a quite special stochastic graph generator, it models an
 * evolving graph. In each time step a single vertex is added to the
 * network and it cites a number of other vertices (as specified by
 * the \p edges_per_step argument). The cited vertices are selected
 * based on the last time they were cited. Time is measured by the
 * addition of vertices and it is binned into \p agebins bins.
 * So if the current time step is \c t and the last citation to a
 * given \c i vertex was made in time step \c t0, then \c
 * (t-t0)/binwidth is calculated where binwidth is \c nodes/agebins+1,
 * in the last expression '/' denotes integer division, so the
 * fraction part is omitted.
 *
 * </para><para>
 * The \p preference argument specifies the preferences for the
 * citation lags, i.e. its first elements contains the attractivity
 * of the very recently cited vertices, etc. The last element is
 * special, it contains the attractivity of the vertices which were
 * never cited. This element should be bigger than zero.
 *
 * </para><para>
 * Note that this function generates networks with multiple edges if
 * \p edges_per_step is bigger than one, call \ref igraph_simplify()
 * on the result to get rid of these edges.
 * \param graph Pointer to an uninitialized graph object, the result
 *     will be stored here.
 * \param node The number of vertices in the network.
 * \param edges_per_node The number of edges to add in each time
 *     step.
 * \param agebins The number of age bins to use.
 * \param preference Pointer to an initialized vector of length
 *     \c agebins+1. This contains the `attractivity' of the various
 *     age bins, the last element is the attractivity of the vertices
 *     which were never cited, and it should be greater than zero.
 *     It is a good idea to have all positive values in this vector.
 *     Preferences cannot be negative.
 * \param directed Logical constant, whether to create directed
 *      networks.
 * \return Error code.
 *
 * \sa \ref igraph_barabasi_aging_game().
 *
 * Time complexity: O(|V|*a+|E|*log|V|), |V| is the number of vertices,
 * |E| is the total number of edges, a is the \p agebins parameter.
 */
int igraph_lastcit_game(igraph_t *graph,
                        igraph_integer_t nodes, igraph_integer_t edges_per_node,
                        igraph_integer_t agebins,
                        const igraph_vector_t *preference,
                        igraph_bool_t directed) {

    long int no_of_nodes = nodes;
    igraph_psumtree_t sumtree;
    igraph_vector_t edges;
    long int i, j, k;
    long int *lastcit;
    long int *index;
    long int binwidth;

    if (agebins != igraph_vector_size(preference) - 1) {
        IGRAPH_ERRORF("The `preference' vector should be of length `agebins' plus one."
                     "Number of agebins is %"IGRAPH_PRId", preference vector is of length %"IGRAPH_PRId"",
                     IGRAPH_EINVAL,
                     agebins, igraph_vector_size(preference));
    }
    if (nodes < 0 ) {
        IGRAPH_ERRORF("Number of nodes should be non-negative, received %"IGRAPH_PRId".",
                     IGRAPH_EINVAL,
                     nodes);
    }
    if (agebins < 1 ) {
        IGRAPH_ERRORF("Number of age bins should be at least 1, received %"IGRAPH_PRId".",
                     IGRAPH_EINVAL,
                     agebins);
    }
    if (VECTOR(*preference)[agebins] <= 0) {
        IGRAPH_ERRORF("The last element of the `preference' vector needs to be positive, but is %g.",
                     IGRAPH_EINVAL,
                     VECTOR(*preference)[agebins]);
    }
    if (igraph_vector_min(preference) < 0) {
        IGRAPH_ERRORF("The preference vector must contain only non-negative values, but found %g.",
                     IGRAPH_EINVAL,
                     igraph_vector_min(preference));
    }

    if (nodes == 0) {
        IGRAPH_CHECK(igraph_empty(graph, nodes, directed));
        return IGRAPH_SUCCESS;
    }

    binwidth = no_of_nodes / agebins + 1;
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    lastcit = IGRAPH_CALLOC(no_of_nodes, long int);
    if (!lastcit) {
        IGRAPH_ERROR("lastcit game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, lastcit);

    index = IGRAPH_CALLOC(no_of_nodes + 1, long int);
    if (!index) {
        IGRAPH_ERROR("lastcit game failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, index);

    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes * edges_per_node));

    /* The first node */
    IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, VECTOR(*preference)[agebins]));
    index[0] = 0;
    index[1] = 0;

    RNG_BEGIN();

    for (i = 1; i < no_of_nodes; i++) {

        /* Add new edges */
        for (j = 0; j < edges_per_node; j++) {
            long int to;
            igraph_real_t sum = igraph_psumtree_sum(&sumtree);
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to);
            lastcit[to] = i + 1;
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, to, VECTOR(*preference)[0]));
        }

        /* Add the node itself */
        IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, VECTOR(*preference)[agebins]));
        index[i + 1] = index[i] + edges_per_node;

        /* Update the preference of some vertices if they got to another bin.
           We need to know the citations of some older vertices, this is in the index. */
        for (k = 1; i - binwidth * k >= 1; k++) {
            long int shnode = i - binwidth * k;
            long int m = index[shnode], n = index[shnode + 1];
            for (j = 2 * m; j < 2 * n; j += 2) {
                long int cnode = (long int) VECTOR(edges)[j + 1];
                if (lastcit[cnode] == shnode + 1) {
                    IGRAPH_CHECK(igraph_psumtree_update(&sumtree, cnode, VECTOR(*preference)[k]));
                }
            }
        }

    }

    RNG_END();

    igraph_psumtree_destroy(&sumtree);
    igraph_free(index);
    igraph_free(lastcit);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_cited_type_game
 * \brief Simulates a citation based on vertex types.
 *
 * Function to create a network based on some vertex categories. This
 * function creates a citation network: in each step a single vertex
 * and \p edges_per_step citing edges are added. Nodes with
 * different categories may have different probabilities to get
 * cited, as given by the \p pref vector.
 *
 * </para><para>
 * Note that this function might generate networks with multiple edges
 * if \p edges_per_step is greater than one. You might want to call
 * \ref igraph_simplify() on the result to remove multiple edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the network.
 * \param types Numeric vector giving the categories of the vertices,
 *     so it should contain \p nodes non-negative integer
 *     numbers. Types are numbered from zero.
 * \param pref The attractivity of the different vertex categories in
 *     a vector. Its length should be the maximum element in \p types
 *     plus one (types are numbered from zero).
 * \param edges_per_step Integer constant, the number of edges to add
 *     in each time step.
 * \param directed Logical constant, whether to create a directed
 *     network.
 * \return Error code.
 *
 * \sa \ref igraph_citing_cited_type_game() for a bit more general
 * game.
 *
 * Time complexity: O((|V|+|E|)log|V|), |V| and |E| are number of
 * vertices and edges, respectively.
 */

int igraph_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
                           const igraph_vector_t *types,
                           const igraph_vector_t *pref,
                           igraph_integer_t edges_per_step,
                           igraph_bool_t directed) {

    igraph_vector_t edges;
    igraph_vector_t cumsum;
    igraph_real_t sum, nnval;
    long int i, j, type;
    long int pref_len = igraph_vector_size(pref);

    if (igraph_vector_size(types) != nodes) {
        IGRAPH_ERRORF("Length of types vector (%ld) must match number of nodes (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, (long) igraph_vector_size(types), nodes);
    }

    if (nodes == 0) {
        igraph_empty(graph, 0, directed);
        return IGRAPH_SUCCESS;
    }

    /* the case of zero-length type vector is caught above, safe to call vector_min here */
    if (igraph_vector_min(types) < 0) {
        IGRAPH_ERRORF("Types should be non-negative, but found %g.",
                      IGRAPH_EINVAL, igraph_vector_min(types));
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    IGRAPH_VECTOR_INIT_FINALLY(&cumsum, 2);
    IGRAPH_CHECK(igraph_vector_reserve(&cumsum, nodes + 1));
    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes * edges_per_step));

    /* first node */
    VECTOR(cumsum)[0] = 0;
    type = (long int) VECTOR(*types)[0];
    if (type >= pref_len) {
        goto err_pref_too_short;
    }
    nnval = VECTOR(*pref)[type];
    if (nnval < 0) {
        goto err_pref_neg;
    }
    sum = VECTOR(cumsum)[1] = nnval;

    RNG_BEGIN();

    for (i = 1; i < nodes; i++) {
        for (j = 0; j < edges_per_step; j++) {
            long int to;
            if (sum > 0) {
                igraph_vector_binsearch(&cumsum, RNG_UNIF(0, sum), &to);
            } else {
                to = i + 1;
            }
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to - 1);
        }
        type = (long int) VECTOR(*types)[i];
        if (type >= pref_len) {
            goto err_pref_too_short;
        }
        nnval = VECTOR(*pref)[type];
        if (nnval < 0) {
            goto err_pref_neg;
        }
        sum += nnval;
        igraph_vector_push_back(&cumsum, sum);
    }

    RNG_END();

    igraph_vector_destroy(&cumsum);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;

err_pref_too_short:
    IGRAPH_ERRORF("Preference vector should have length at least %ld with the given types.", IGRAPH_EINVAL,
                  (long) igraph_vector_max(types) + 1);

err_pref_neg:
    IGRAPH_ERRORF("Preferences should be non-negative, but found %g.", IGRAPH_EINVAL,
                  igraph_vector_min(pref));
}

static void igraph_i_citing_cited_type_game_free(igraph_i_citing_cited_type_game_struct_t *s) {
    long int i;
    if (!s->sumtrees) {
        return;
    }
    for (i = 0; i < s->no; i++) {
        igraph_psumtree_destroy(&s->sumtrees[i]);
    }
    igraph_free(s->sumtrees);
}

/**
 * \function igraph_citing_cited_type_game
 * \brief Simulates a citation network based on vertex types.
 *
 * This game is similar to \ref igraph_cited_type_game() but here the
 * category of the citing vertex is also considered.
 *
 * </para><para>
 * An evolving citation network is modeled here, a single vertex and
 * its \p edges_per_step citation are added in each time step. The
 * odds the a given vertex is cited by the new vertex depends on the
 * category of both the citing and the cited vertex and is given in
 * the \p pref matrix. The categories of the citing vertex correspond
 * to the rows, the categories of the cited vertex to the columns of
 * this matrix. I.e. the element in row \c i and column \c j gives the
 * probability that a \c j vertex is cited, if the category of the
 * citing vertex is \c i.
 *
 * </para><para>
 * Note that this function might generate networks with multiple edges
 * if \p edges_per_step is greater than one. You might want to call
 * \ref igraph_simplify() on the result to remove multiple edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the network.
 * \param types A numeric matrix of length \p nodes, containing the
 *    categories of the vertices. The categories are numbered from
 *    zero.
 * \param pref The preference matrix, a square matrix is required,
 *     both the number of rows and columns should be the maximum
 *     element in \p types plus one (types are numbered from zero).
 * \param directed Logical constant, whether to create a directed
 *     network.
 * \return Error code.
 *
 * Time complexity: O((|V|+|E|)log|V|), |V| and |E| are number of
 * vertices and edges, respectively.
 */

int igraph_citing_cited_type_game(igraph_t *graph, igraph_integer_t nodes,
                                  const igraph_vector_t *types,
                                  const igraph_matrix_t *pref,
                                  igraph_integer_t edges_per_step,
                                  igraph_bool_t directed) {

    igraph_vector_t edges;
    igraph_i_citing_cited_type_game_struct_t str = { 0, NULL };
    igraph_psumtree_t *sumtrees;
    igraph_vector_t sums;
    long int no_of_types;
    long int i, j;

    if (igraph_vector_size(types) != nodes) {
        IGRAPH_ERRORF("Length of types vector (%ld) not equal to number"
                      " of nodes (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(types), nodes);
    }

    /* avoid calling vector_max on empty vector */
    no_of_types = nodes == 0 ? 0 : igraph_vector_max(types) + 1;

    if (igraph_matrix_ncol(pref) != no_of_types) {
        IGRAPH_ERRORF("Number of preference matrix colums (%ld) not "
                      "equal to number of types (%g).",
                      IGRAPH_EINVAL,
                      igraph_matrix_ncol(pref),
                      no_of_types);
    }
    if (igraph_matrix_nrow(pref) != no_of_types) {
        IGRAPH_ERRORF("Number of preference matrix rows (%ld) not "
                      "equal to number of types (%g).",
                      IGRAPH_EINVAL,
                      igraph_matrix_nrow(pref),
                      no_of_types);
    }

    /* return an empty graph if nodes is zero */
    if (nodes == 0) {
        return igraph_empty(graph, 0, directed);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

    str.sumtrees = sumtrees = IGRAPH_CALLOC(no_of_types, igraph_psumtree_t);
    if (!sumtrees) {
        IGRAPH_ERROR("Citing-cited type game failed.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_citing_cited_type_game_free, &str);

    for (i = 0; i < no_of_types; i++) {
        IGRAPH_CHECK(igraph_psumtree_init(&sumtrees[i], nodes));
        str.no++;
    }
    IGRAPH_VECTOR_INIT_FINALLY(&sums, no_of_types);

    IGRAPH_CHECK(igraph_vector_reserve(&edges, nodes * edges_per_step));

    /* First node */
    for (i = 0; i < no_of_types; i++) {
        long int type = (long int) VECTOR(*types)[0];
        if ( MATRIX(*pref, i, type) < 0) {
            IGRAPH_ERRORF("Preference matrix contains negative entry: %g.", IGRAPH_EINVAL, MATRIX(*pref, i, type));
        }
        IGRAPH_CHECK(igraph_psumtree_update(&sumtrees[i], 0, MATRIX(*pref, i, type)));
        VECTOR(sums)[i] = MATRIX(*pref, i, type);
    }

    RNG_BEGIN();

    for (i = 1; i < nodes; i++) {
        long int type = (long int) VECTOR(*types)[i];
        igraph_real_t sum = VECTOR(sums)[type];
        for (j = 0; j < edges_per_step; j++) {
            long int to;
            igraph_psumtree_search(&sumtrees[type], &to, RNG_UNIF(0, sum));
            igraph_vector_push_back(&edges, i);
            igraph_vector_push_back(&edges, to);
        }

        /* add i */
        for (j = 0; j < no_of_types; j++) {
            if ( MATRIX(*pref, j, type) < 0) {
                IGRAPH_ERRORF("Preference matrix contains negative entry: %g.", IGRAPH_EINVAL, MATRIX(*pref, j, type));
            }
            IGRAPH_CHECK(igraph_psumtree_update(&sumtrees[j], i, MATRIX(*pref, j,  type)));
            VECTOR(sums)[j] += MATRIX(*pref, j, type);
        }
    }

    RNG_END();

    igraph_i_citing_cited_type_game_free(&str);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_create(graph, &edges, nodes, directed);
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&sums);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}
