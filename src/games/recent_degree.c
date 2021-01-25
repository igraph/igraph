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
#include "igraph_dqueue.h"
#include "igraph_psumtree.h"
#include "igraph_random.h"

/**
 * \function igraph_recent_degree_game
 * \brief Stochastic graph generator based on the number of incident edges a node has gained recently
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param n The number of vertices in the graph, this is the same as
 *        the number of time steps.
 * \param power The exponent, the probability that a node gains a
 *        new edge is proportional to the number of edges it has
 *        gained recently (in the last \p window time steps) to \p
 *        power.
 * \param window Integer constant, the size of the time window to use
 *        to count the number of recent edges.
 * \param m Integer constant, the number of edges to add per time
 *        step if the \p outseq parameter is a null pointer or a
 *        zero-length vector.
 * \param outseq The number of edges to add in each time step. This
 *        argument is ignored if it is a null pointer or a zero length
 *        vector, is this case the constant \p m parameter is used.
 * \param outpref Logical constant, if true the edges originated by a
 *        vertex also count as recent incident edges. It is false in
 *        most cases.
 * \param zero_appeal Constant giving the attractiveness of the
 *        vertices which haven't gained any edge recently.
 * \param directed Logical constant, whether to generate a directed
 *        graph.
 * \return Error code.
 *
 * Time complexity: O(|V|*log(|V|)+|E|), |V| is the number of
 * vertices, |E| is the number of edges in the graph.
 *
 */
int igraph_recent_degree_game(igraph_t *graph, igraph_integer_t n,
                              igraph_real_t power,
                              igraph_integer_t window,
                              igraph_integer_t m,
                              const igraph_vector_t *outseq,
                              igraph_bool_t outpref,
                              igraph_real_t zero_appeal,
                              igraph_bool_t directed) {

    long int no_of_nodes = n;
    long int no_of_neighbors = m;
    long int no_of_edges;
    igraph_vector_t edges;
    long int i, j;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;
    long int time_window = window;
    igraph_dqueue_t history;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != n) {
        IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
    }

    if (outseq == 0 || igraph_vector_size(outseq) == 0) {
        no_of_neighbors = m;
        no_of_edges = (no_of_nodes - 1) * no_of_neighbors;
    } else {
        no_of_edges = 0;
        for (i = 1; i < igraph_vector_size(outseq); i++) {
            no_of_edges += VECTOR(*outseq)[i];
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_dqueue_init(&history,
                                    time_window * (no_of_neighbors + 1) + 10));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &history);

    RNG_BEGIN();

    /* first node */
    IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, zero_appeal));
    igraph_dqueue_push(&history, -1);

    /* and the rest */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        long int to;
        if (outseq != 0 && igraph_vector_size(outseq) != 0) {
            no_of_neighbors = (long int) VECTOR(*outseq)[i];
        }

        if (i >= time_window) {
            while ((j = (long int) igraph_dqueue_pop(&history)) != -1) {
                VECTOR(degree)[j] -= 1;
                IGRAPH_CHECK(igraph_psumtree_update(&sumtree, j, pow(VECTOR(degree)[j], power) + zero_appeal));
            }
        }

        sum = igraph_psumtree_sum(&sumtree);
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
            igraph_dqueue_push(&history, to);
        }
        igraph_dqueue_push(&history, -1);

        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            long int nn = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, nn, pow(VECTOR(degree)[nn], power) + zero_appeal));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, pow(VECTOR(degree)[i], power) + zero_appeal));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, zero_appeal));
        }
    }

    RNG_END();

    igraph_dqueue_destroy(&history);
    igraph_psumtree_destroy(&sumtree);
    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_recent_degree_aging_game
 * \brief Preferential attachment based on the number of edges gained recently, with aging of vertices
 *
 * </para><para>
 * This game is very similar to \ref igraph_barabasi_aging_game(),
 * except that instead of the total number of incident edges the
 * number of edges gained in the last \p time_window time steps are
 * counted.
 *
 * </para><para>The degree dependent part of the attractiveness is
 * given by k to the power of \p pa_exp plus \p zero_appeal; the age
 * dependent part is l to the power to \p aging_exp.
 * k is the number of edges gained in the last \p time_window time
 * steps, l is the age of the vertex.
 * \param graph Pointer to an uninitialized graph object.
 * \param nodes The number of vertices in the graph.
 * \param m The number of edges to add in each time step. If the \p
 *        outseq argument is not a null vector or a zero-length vector
 *        then it is ignored.
 * \param outseq Vector giving the number of edges to add in each time
 *        step. If it is a null pointer or a zero-length vector then
 *        it is ignored and the \p m argument is used.
 * \param outpref Logical constant, if true the edges initiated by a
 *        vertex are also counted. Normally it is false.
 * \param pa_exp The exponent for the preferential attachment.
 * \param aging_exp The exponent for the aging, normally it is
 *        negative: old vertices gain edges with less probability.
 * \param aging_bin Integer constant, gives the scale of the aging.
 *        The age of the vertices is incremented by one after every \p
 *        aging_bin vertex added.
 * \param time_window The time window to use to count the number of
 *        incident edges for the vertices.
 * \param zero_appeal The degree dependent part of the attractiveness
 *        for zero degree vertices.
 * \param directed Logical constant, whether to create a directed
 *        graph.
 * \return Error code.
 *
 * Time complexity: O((|V|+|V|/aging_bin)*log(|V|)+|E|). |V| is the number
 * of vertices, |E| the number of edges.
 */
int igraph_recent_degree_aging_game(igraph_t *graph,
                                    igraph_integer_t nodes,
                                    igraph_integer_t m,
                                    const igraph_vector_t *outseq,
                                    igraph_bool_t outpref,
                                    igraph_real_t pa_exp,
                                    igraph_real_t aging_exp,
                                    igraph_integer_t aging_bin,
                                    igraph_integer_t time_window,
                                    igraph_real_t zero_appeal,
                                    igraph_bool_t directed) {

    long int no_of_nodes = nodes;
    long int no_of_neighbors = m;
    long int binwidth = nodes / aging_bin + 1;
    long int no_of_edges;
    igraph_vector_t edges;
    long int i, j, k;
    igraph_psumtree_t sumtree;
    long int edgeptr = 0;
    igraph_vector_t degree;
    igraph_dqueue_t history;

    if (no_of_nodes < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (outseq != 0 && igraph_vector_size(outseq) != 0 && igraph_vector_size(outseq) != no_of_nodes) {
        IGRAPH_ERROR("Invalid out degree sequence length", IGRAPH_EINVAL);
    }
    if ( (outseq == 0 || igraph_vector_size(outseq) == 0) && m < 0) {
        IGRAPH_ERROR("Invalid out degree", IGRAPH_EINVAL);
    }
    if (aging_bin <= 0) {
        IGRAPH_ERROR("Invalid aging bin", IGRAPH_EINVAL);
    }

    if (outseq == 0 || igraph_vector_size(outseq) == 0) {
        no_of_neighbors = m;
        no_of_edges = (no_of_nodes - 1) * no_of_neighbors;
    } else {
        no_of_edges = 0;
        for (i = 1; i < igraph_vector_size(outseq); i++) {
            no_of_edges += VECTOR(*outseq)[i];
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_psumtree_init(&sumtree, no_of_nodes));
    IGRAPH_FINALLY(igraph_psumtree_destroy, &sumtree);
    IGRAPH_VECTOR_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_dqueue_init(&history,
                                    time_window * (no_of_neighbors + 1) + 10));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &history);

    RNG_BEGIN();

    /* first node */
    IGRAPH_CHECK(igraph_psumtree_update(&sumtree, 0, zero_appeal));
    igraph_dqueue_push(&history, -1);

    /* and the rest */
    for (i = 1; i < no_of_nodes; i++) {
        igraph_real_t sum;
        long int to;
        if (outseq != 0 && igraph_vector_size(outseq) != 0) {
            no_of_neighbors = (long int) VECTOR(*outseq)[i];
        }

        if (i >= time_window) {
            while ((j = (long int) igraph_dqueue_pop(&history)) != -1) {
                long int age = (i - j) / binwidth;
                VECTOR(degree)[j] -= 1;
                IGRAPH_CHECK(igraph_psumtree_update(
                    &sumtree, j,
                    (pow(VECTOR(degree)[j], pa_exp) + zero_appeal) * pow(age + 1, aging_exp)
                ));
            }
        }

        sum = igraph_psumtree_sum(&sumtree);
        for (j = 0; j < no_of_neighbors; j++) {
            igraph_psumtree_search(&sumtree, &to, RNG_UNIF(0, sum));
            VECTOR(degree)[to]++;
            VECTOR(edges)[edgeptr++] = i;
            VECTOR(edges)[edgeptr++] = to;
            igraph_dqueue_push(&history, to);
        }
        igraph_dqueue_push(&history, -1);

        /* update probabilities */
        for (j = 0; j < no_of_neighbors; j++) {
            long int n = (long int) VECTOR(edges)[edgeptr - 2 * j - 1];
            long int age = (i - n) / binwidth;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, n,
                (pow(VECTOR(degree)[n], pa_exp) + zero_appeal) * pow(age + 1, aging_exp)
            ));
        }
        if (outpref) {
            VECTOR(degree)[i] += no_of_neighbors;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, i,
                pow(VECTOR(degree)[i], pa_exp) + zero_appeal
            ));
        } else {
            IGRAPH_CHECK(igraph_psumtree_update(&sumtree, i, zero_appeal));
        }

        /* aging */
        for (k = 1; i - binwidth * k + 1 >= 1; k++) {
            long int shnode = i - binwidth * k;
            long int deg = (long int) VECTOR(degree)[shnode];
            long int age = (i - shnode) / binwidth;
            IGRAPH_CHECK(igraph_psumtree_update(
                &sumtree, shnode,
                (pow(deg, pa_exp) + zero_appeal) * pow(age + 2, aging_exp)
            ));
        }
    }

    RNG_END();

    igraph_dqueue_destroy(&history);
    igraph_vector_destroy(&degree);
    igraph_psumtree_destroy(&sumtree);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
