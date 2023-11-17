/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2020 The igraph development team

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

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "core/interruption.h"

#include <limits.h>

/**
 * \section about_leading_eigenvector_methods
 *
 * <para>
 * The function documented in these section implements the
 * <quote>leading eigenvector</quote> method developed by Mark Newman and
 * published in MEJ Newman: Finding community structure using the
 * eigenvectors of matrices, Phys Rev E 74:036104 (2006).</para>
 *
 * <para>
 * The heart of the method is the definition of the modularity matrix,
 * B, which is B=A-P, A being the adjacency matrix of the (undirected)
 * network, and P contains the probability that certain edges are
 * present according to the <quote>configuration model</quote> In
 * other words, a Pij element of P is the probability that there is an
 * edge between vertices i and j in a random network in which the
 * degrees of all vertices are the same as in the input graph.</para>
 *
 * <para>
 * The leading eigenvector method works by calculating the eigenvector
 * of the modularity matrix for the largest positive eigenvalue and
 * then separating vertices into two community based on the sign of
 * the corresponding element in the eigenvector. If all elements in
 * the eigenvector are of the same sign that means that the network
 * has no underlying community structure.
 * Check Newman's paper to understand why this is a good method for
 * detecting community structure. </para>
 *
 * <para>
 * The leading eigenvector community structure detection method is
 * implemented in \ref igraph_community_leading_eigenvector(). After
 * the initial split, the following splits are done in a way to
 * optimize modularity regarding to the original network. Note that
 * any further refinement, for example using Kernighan-Lin, as
 * proposed in Section V.A of Newman (2006), is not implemented here.
 * </para>
 *
 * <para>
 * \example examples/simple/igraph_community_leading_eigenvector.c
 * </para>
 */

typedef struct igraph_i_community_leading_eigenvector_data_t {
    igraph_vector_int_t *idx;
    igraph_vector_int_t *idx2;
    igraph_adjlist_t *adjlist;
    igraph_inclist_t *inclist;
    igraph_vector_t *tmp;
    igraph_integer_t no_of_edges;
    igraph_vector_int_t *mymembership;
    igraph_integer_t comm;
    const igraph_vector_t *weights;
    const igraph_t *graph;
    igraph_vector_t *strength;
    igraph_real_t sumweights;
} igraph_i_community_leading_eigenvector_data_t;

static igraph_error_t igraph_i_community_leading_eigenvector(
        igraph_real_t *to,
        const igraph_real_t *from,
        int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    igraph_integer_t size = n;
    igraph_vector_int_t *idx = data->idx;
    igraph_vector_int_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_real_t ktx, ktx2;
    igraph_integer_t no_of_edges = data->no_of_edges;
    igraph_vector_int_t *mymembership = data->mymembership;
    igraph_integer_t comm = data->comm;

    /* Ax */
    for (igraph_integer_t j = 0; j < size; j++) {
        igraph_integer_t oldid = VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_integer_t nlen = igraph_vector_int_size(neis);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (igraph_integer_t k = 0; k < nlen; k++) {
            igraph_integer_t nei = VECTOR(*neis)[k];
            igraph_integer_t neimemb = VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                to[j] += from[ VECTOR(*idx2)[nei] ];
                VECTOR(*tmp)[j] += 1;
            }
        }
    }

    /* Now calculate k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (igraph_integer_t j = 0; j < size; j++) {
        igraph_integer_t oldid = VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_integer_t degree = igraph_vector_int_size(neis);
        ktx += from[j] * degree;
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges / 2.0;
    ktx2 = ktx2 / no_of_edges / 2.0;

    /* Now calculate Bx */
    for (igraph_integer_t j = 0; j < size; j++) {
        igraph_integer_t oldid = VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree = igraph_vector_int_size(neis);
        to[j] = to[j] - ktx * degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * degree;
    }

    /* -d_ij summa l in G B_il */
    for (igraph_integer_t j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_community_leading_eigenvector_weighted(
        igraph_real_t *to,
        const igraph_real_t *from,
        int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    igraph_integer_t size = n;
    igraph_vector_int_t *idx = data->idx;
    igraph_vector_int_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_inclist_t *inclist = data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_int_t *mymembership = data->mymembership;
    igraph_integer_t comm = data->comm;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *strength = data->strength;
    igraph_real_t sw = data->sumweights;

    /* Ax */
    for (igraph_integer_t j = 0; j < size; j++) {
        igraph_integer_t oldid = VECTOR(*idx)[j];
        igraph_vector_int_t *inc = igraph_inclist_get(inclist, oldid);
        igraph_integer_t nlen = igraph_vector_int_size(inc);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (igraph_integer_t k = 0; k < nlen; k++) {
            igraph_integer_t edge = VECTOR(*inc)[k];
            igraph_real_t w = VECTOR(*weights)[edge];
            igraph_integer_t nei = IGRAPH_OTHER(graph, edge, oldid);
            igraph_integer_t neimemb = VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                to[j] += from[ VECTOR(*idx2)[nei] ] * w;
                VECTOR(*tmp)[j] += w;
            }
        }
    }

    /* k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (igraph_integer_t j = 0; j < size; j++) {
        igraph_integer_t oldid = VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        ktx += from[j] * str;
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;

    /* Bx */
    for (igraph_integer_t j = 0; j < size; j++) {
        igraph_integer_t oldid = VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * str;
    }

    /* -d_ij summa l in G B_il */
    for (igraph_integer_t j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return IGRAPH_SUCCESS;
}

static void igraph_i_error_handler_none(const char *reason, const char *file,
                                        int line, igraph_error_t igraph_errno) {
    IGRAPH_UNUSED(reason);
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(igraph_errno);
    /* do nothing */
}


/**
 * \ingroup communities
 * \function igraph_community_leading_eigenvector
 * \brief Leading eigenvector community finding (proper version).
 *
 * Newman's leading eigenvector method for detecting community
 * structure. This is the proper implementation of the recursive,
 * divisive algorithm: each split is done by maximizing the modularity
 * regarding the original network, see MEJ Newman: Finding community
 * structure in networks using the eigenvectors of matrices,
 * Phys Rev E 74:036104 (2006).
 * https://doi.org/10.1103/PhysRevE.74.036104
 *
 * \param graph The input graph. Edge directions will be ignored.
 * \param weights The weights of the edges, or a null pointer for
 *    unweighted graphs.
 * \param merges The result of the algorithm, a matrix containing the
 *    information about the splits performed. The matrix is built in
 *    the opposite way however, it is like the result of an
 *    agglomerative algorithm. Unlike with most other hierarchicaly
 *    community detection functions in igraph, the integers in this matrix
 *    represent community indices, not vertex indices. If at the end of
 *    the algorithm (after \p steps steps was done) there are <quote>p</quote>
 *    communities, then these are numbered from zero to <quote>p-1</quote>.
 *    The first line of the matrix contains the first <quote>merge</quote>
 *    (which is in reality the last split) of two communities into
 *    community <quote>p</quote>, the merge in the second line forms
 *    community <quote>p+1</quote>, etc. The matrix should be
 *    initialized before calling and will be resized as needed.
 *    This argument is ignored if it is \c NULL.
 * \param membership The membership of the vertices after all the
 *    splits were performed will be stored here. The vector must be
 *    initialized  before calling and will be resized as needed.
 *    This argument is ignored if it is \c NULL. This argument can
 *    also be used to supply a starting configuration for the community
 *    finding, in the format of a membership vector. In this case the
 *    \p start argument must be set to 1.
 * \param steps The maximum number of steps to perform. It might
 *    happen that some component (or the whole network) has no
 *    underlying community structure and no further steps can be
 *    done. If you want as many steps as possible then supply the
 *    number of vertices in the network here.
 * \param options The options for ARPACK. Supply \c NULL here to use the
 *    defaults. \c n is always overwritten. \c ncv is set to at least 4.
 * \param modularity If not a null pointer, then it must be a pointer
 *    to a real number and the modularity score of the final division
 *    is stored here.
 * \param start Boolean, whether to use the community structure given
 *    in the \p membership argument as a starting point.
 * \param eigenvalues Pointer to an initialized vector or a null
 *    pointer. If not a null pointer, then the eigenvalues calculated
 *    along the community structure detection are stored here. The
 *    non-positive eigenvalues, that do not result a split, are stored
 *    as well.
 * \param eigenvectors If not a null pointer, then the eigenvectors
 *    that are calculated in each step of the algorithm are stored here,
 *    in a list of vectors. Each eigenvector is stored in an
 *    \ref igraph_vector_t object.
 * \param history Pointer to an initialized vector or a null pointer.
 *    If not a null pointer, then a trace of the algorithm is stored
 *    here, encoded numerically. The various operations:
 *    \clist
 *    \cli IGRAPH_LEVC_HIST_START_FULL
 *      Start the algorithm from an initial state where each connected
 *      component is a separate community.
 *    \cli IGRAPH_LEVC_HIST_START_GIVEN
 *      Start the algorithm from a given community structure. The next
 *      value in the vector contains the initial number of
 *      communities.
 *    \cli IGRAPH_LEVC_HIST_SPLIT
 *      Split a community into two communities. The id of the splitted
 *      community is given in the next element of the history vector.
 *      The id of the first new community is the same as the id of the
 *      splitted community. The id of the second community equals to
 *      the number of communities before the split.
 *    \cli IGRAPH_LEVC_HIST_FAILED
 *      Tried to split a community, but it was not worth it, as it
 *      does not result in a bigger modularity value. The id of the
 *      community is given in the next element of the vector.
 *    \endclist
 * \param callback A null pointer or a function of type \ref
 *    igraph_community_leading_eigenvector_callback_t. If given, this
 *    callback function is called after each eigenvector/eigenvalue
 *    calculation. If the callback returns \c IGRAPH_STOP, then the
 *    community finding algorithm stops. If it returns \c IGRAPH_SUCCESS,
 *    the algorithm continues normally. Any other return value is considered
 *    an igraph error code and will terminete the algorithm with the same
 *    error code. See the arguments passed to the callback at the documentation
 *    of \ref igraph_community_leading_eigenvector_callback_t.
 * \param callback_extra Extra argument to pass to the callback
 *    function.
 * \return Error code.
 *
 * \sa \ref igraph_community_walktrap() and \ref
 * igraph_community_spinglass() for other community structure
 * detection methods.
 *
 * Time complexity: O(|E|+|V|^2*steps), |V| is the number of vertices,
 * |E| the number of edges, <quote>steps</quote> the number of splits
 * performed.
 */
igraph_error_t igraph_community_leading_eigenvector(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_matrix_int_t *merges,
        igraph_vector_int_t *membership,
        igraph_integer_t steps,
        igraph_arpack_options_t *options,
        igraph_real_t *modularity,
        igraph_bool_t start,
        igraph_vector_t *eigenvalues,
        igraph_vector_list_t *eigenvectors,
        igraph_vector_t *history,
        igraph_community_leading_eigenvector_callback_t *callback,
        void *callback_extra) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_dqueue_int_t tosplit;
    igraph_vector_int_t idx, idx2;
    igraph_vector_t mymerges;
    igraph_vector_t strength, tmp;
    igraph_vector_t start_vec;
    igraph_integer_t staken = 0;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    igraph_integer_t i, j, k, l;
    igraph_integer_t communities;
    igraph_vector_int_t vmembership, *mymembership = membership;
    igraph_i_community_leading_eigenvector_data_t extra;
    igraph_arpack_storage_t storage;
    igraph_real_t mod = 0;
    igraph_arpack_function_t *arpcb1 =
            weights ? igraph_i_community_leading_eigenvector_weighted :
                      igraph_i_community_leading_eigenvector;
    igraph_real_t sumweights = 0.0;

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph too large for ARPACK.", IGRAPH_EOVERFLOW);
    }

    if (weights && no_of_edges != igraph_vector_size(weights)) {
        IGRAPH_ERROR("Weight vector length does not match number of edges.", IGRAPH_EINVAL);
    }

    if (start && !membership) {
        IGRAPH_ERROR("Cannot start from given configuration if memberships missing.", IGRAPH_EINVAL);
    }

    if (start && membership &&
        igraph_vector_int_size(membership) != no_of_nodes) {
        IGRAPH_ERROR("Supplied memberhsip vector length does not match number of vertices.",
                     IGRAPH_EINVAL);
    }

    if (start && membership && igraph_vector_int_max(membership) >= no_of_nodes) {
        IGRAPH_WARNING("Too many communities in membership start vector.");
    }

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Directed graph supplied, edge directions will be ignored.");
    }

    if (steps < 0 || steps > no_of_nodes - 1) {
        steps = no_of_nodes > 0 ? no_of_nodes - 1 : 0;
    }

    if (!membership) {
        mymembership = &vmembership;
        IGRAPH_VECTOR_INT_INIT_FINALLY(mymembership, 0);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&mymerges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&mymerges, steps * 2));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&idx, 0);
    if (eigenvalues)  {
        igraph_vector_clear(eigenvalues);
    }
    if (eigenvectors) {
        igraph_vector_list_clear(eigenvectors);
    }

    if (!start) {
        /* Calculate the weakly connected components in the graph and use them as
         * an initial split */
        IGRAPH_CHECK(igraph_connected_components(graph, mymembership, &idx, 0, IGRAPH_WEAK));
        communities = igraph_vector_int_size(&idx);
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history,
                                                 IGRAPH_LEVC_HIST_START_FULL));
        }
    } else {
        /* Just create the idx vector for the given membership vector */
        communities = igraph_vector_int_max(mymembership) + 1;
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history,
                                                 IGRAPH_LEVC_HIST_START_GIVEN));
            IGRAPH_CHECK(igraph_vector_push_back(history, communities));
        }
        IGRAPH_CHECK(igraph_vector_int_resize(&idx, communities));
        igraph_vector_int_null(&idx);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_integer_t t = VECTOR(*mymembership)[i];
            VECTOR(idx)[t] += 1;
        }
    }

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&tosplit, 100);
    for (i = 0; i < communities; i++) {
        if (VECTOR(idx)[i] > 2) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&tosplit, i));
        }
    }
    for (i = 1; i < communities; i++) {
        /* Record merge */
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, i - 1));
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, i));
        if (eigenvalues) {
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, IGRAPH_NAN));
        }
        if (eigenvectors) {
            /* There are no eigenvectors associated to these steps because the
             * splits were given by the user (or by the components of the graph)
             * so we push empty vectors */
            IGRAPH_CHECK(igraph_vector_list_push_back_new(eigenvectors, NULL));
        }
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, i - 1));
        }
    }
    staken = communities - 1;

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_int_resize(&idx, no_of_nodes));
    igraph_vector_int_null(&idx);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&idx2, no_of_nodes);
    if (!weights) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    } else {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
        IGRAPH_VECTOR_INIT_FINALLY(&strength, no_of_nodes);
        IGRAPH_CHECK(igraph_strength(graph, &strength, igraph_vss_all(),
                                     IGRAPH_ALL, IGRAPH_LOOPS, weights));
        sumweights = igraph_vector_sum(weights);
    }

    if (options == NULL) {
        options = igraph_arpack_options_get_default();
    }

    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->which[0] = 'L'; options->which[1] = 'A';

    /* Memory for ARPACK */
    /* We are allocating memory for 20 eigenvectors since options->ncv won't be
     * larger than 20 when using automatic mode in igraph_arpack_rssolve */
    IGRAPH_CHECK(igraph_arpack_storage_init(&storage, (int) no_of_nodes, 20,
                                            (int) no_of_nodes, 1));
    IGRAPH_FINALLY(igraph_arpack_storage_destroy, &storage);
    extra.idx = &idx;
    extra.idx2 = &idx2;
    extra.tmp = &tmp;
    extra.adjlist = &adjlist;
    extra.inclist = &inclist;
    extra.weights = weights;
    extra.sumweights = sumweights;
    extra.graph = graph;
    extra.strength = &strength;
    extra.no_of_edges = no_of_edges;
    extra.mymembership = mymembership;

    while (!igraph_dqueue_int_empty(&tosplit) && staken < steps) {
        igraph_integer_t comm = igraph_dqueue_int_pop_back(&tosplit);
        /* depth first search */
        igraph_integer_t size = 0;

        IGRAPH_ALLOW_INTERRUPTION();

        for (i = 0; i < no_of_nodes; i++) {
            if (VECTOR(*mymembership)[i] == comm) {
                VECTOR(idx)[size] = i;
                VECTOR(idx2)[i] = size++;
            }
        }

        staken++;
        if (size <= 2) {
            continue;
        }

        options->n = (int) size;
        options->info = 0;
        options->nev = 1;
        options->ldv = 0;
        options->ncv = 0; /* 0 means "automatic" in igraph_arpack_rssolve */
        options->nconv = 0;
        options->lworkl = 0; /* we surely have enough space */
        extra.comm = comm;

        /* Use a random start vector, but don't let ARPACK generate the
         * start vector -- we want to use our own RNG. Also, we want to generate
         * values close to +1 and -1 as this is what the eigenvector should
         * look like if there _is_ some kind of a community structure at this
         * step to discover. Experiments showed that shuffling a vector
         * containing equal number of slightly perturbed +/-1 values yields
         * convergence in most cases. */
        options->start = 1;
        options->mxiter = options->mxiter > 10000 ? options->mxiter : 10000;  /* use more iterations, we've had convergence problems with 3000 */
        RNG_BEGIN();
        for (i = 0; i < options->n; i++) {
            storage.resid[i] = (i % 2 ? 1 : -1) + RNG_UNIF(-0.1, 0.1);
        }
        RNG_END();
        igraph_vector_view(&start_vec, storage.resid, options->n);
        IGRAPH_CHECK(igraph_vector_shuffle(&start_vec));

        {
            igraph_error_t retval;
            igraph_error_handler_t *errh =
                    igraph_set_error_handler(igraph_i_error_handler_none);
            retval = igraph_arpack_rssolve(arpcb1, &extra, options, &storage, /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            if (retval != IGRAPH_SUCCESS && retval != IGRAPH_ARPACK_MAXIT && retval != IGRAPH_ARPACK_NOSHIFT) {
                IGRAPH_ERROR("ARPACK call failed", retval);
            }
        }

        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }

        /* Ok, we have the leading eigenvector of the modularity matrix */

        /* ---------------------------------------------------------------*/
        /* To avoid numeric errors */
        if (fabs(storage.d[0]) < 1e-8) {
            storage.d[0] = 0;
        }

        /* We replace very small (in absolute value) elements of the
           leading eigenvector with zero, to get the same result,
           consistently.*/
        for (i = 0; i < size; i++) {
            if (fabs(storage.v[i]) < 1e-8) {
                storage.v[i] = 0;
            }
        }

        /* Just to have the always the same result, we multiply by -1
           if the first (nonzero) element is not positive. */
        for (i = 0; i < size; i++) {
            if (storage.v[i] != 0) {
                break;
            }
        }
        if (i < size && storage.v[i] < 0) {
            for (i = 0; i < size; i++) {
                storage.v[i] = - storage.v[i];
            }
        }
        /* ---------------------------------------------------------------*/

        if (callback) {
            igraph_vector_t vv;
            igraph_error_t ret;

            igraph_vector_view(&vv, storage.v, size);
            IGRAPH_CHECK_CALLBACK(
                        callback(
                            mymembership, comm, storage.d[0], &vv, arpcb1,
                            &extra, callback_extra
                        ), &ret
                    );

            if (ret == IGRAPH_STOP) {
                break;
            }
        }

        if (eigenvalues) {
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, storage.d[0]));
        }

        if (eigenvectors) {
            igraph_vector_t *v;
            /* TODO: this would be faster if we had an igraph_vector_list_push_back_new_with_size_hint */
            IGRAPH_CHECK(igraph_vector_list_push_back_new(eigenvectors, &v));
            IGRAPH_CHECK(igraph_vector_resize(v, size));
            for (i = 0; i < size; i++) {
                VECTOR(*v)[i] = storage.v[i];
            }
        }

        if (storage.d[0] <= 0) {
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        /* Count the number of vertices in each community after the split */
        l = 0;
        for (j = 0; j < size; j++) {
            if (storage.v[j] < 0) {
                storage.v[j] = -1;
                l++;
            } else {
                storage.v[j] = 1;
            }
        }
        if (l == 0 || l == size) {
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        /* Check that Q increases with our choice of split */
        arpcb1(storage.v + size, storage.v, (int) size, &extra);
        mod = 0;
        for (i = 0; i < size; i++) {
            mod += storage.v[size + i] * storage.v[i];
        }
        if (mod <= 1e-8) {
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        communities++;

        /* Rewrite the mymembership vector */
        for (j = 0; j < size; j++) {
            if (storage.v[j] < 0) {
                igraph_integer_t oldid = VECTOR(idx)[j];
                VECTOR(*mymembership)[oldid] = communities - 1;
            }
        }

        /* Record merge */
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, comm));
        IGRAPH_CHECK(igraph_vector_push_back(&mymerges, communities - 1));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, comm));
        }

        /* Store the resulting communities in the queue if needed */
        if (l > 1) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&tosplit, communities - 1));
        }
        if (size - l > 1) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&tosplit, comm));
        }

    }

    igraph_arpack_storage_destroy(&storage);
    IGRAPH_FINALLY_CLEAN(1);
    if (!weights) {
        igraph_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_inclist_destroy(&inclist);
        igraph_vector_destroy(&strength);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_dqueue_int_destroy(&tosplit);
    igraph_vector_destroy(&tmp);
    igraph_vector_int_destroy(&idx2);
    IGRAPH_FINALLY_CLEAN(3);

    /* reform the mymerges vector */
    if (merges) {
        igraph_vector_int_null(&idx);
        l = igraph_vector_size(&mymerges);
        k = communities;
        j = 0;
        IGRAPH_CHECK(igraph_matrix_int_resize(merges, l / 2, 2));
        for (i = l; i > 0; i -= 2) {
            igraph_integer_t from = VECTOR(mymerges)[i - 1];
            igraph_integer_t to = VECTOR(mymerges)[i - 2];
            MATRIX(*merges, j, 0) = VECTOR(mymerges)[i - 2];
            MATRIX(*merges, j, 1) = VECTOR(mymerges)[i - 1];
            if (VECTOR(idx)[from] != 0) {
                MATRIX(*merges, j, 1) = VECTOR(idx)[from] - 1;
            }
            if (VECTOR(idx)[to] != 0) {
                MATRIX(*merges, j, 0) = VECTOR(idx)[to] - 1;
            }
            VECTOR(idx)[to] = ++k;
            j++;
        }
    }

    igraph_vector_int_destroy(&idx);
    igraph_vector_destroy(&mymerges);
    IGRAPH_FINALLY_CLEAN(2);

    if (modularity) {
        IGRAPH_CHECK(igraph_modularity(graph, mymembership, weights,
                                       /* resolution */ 1,
                                       IGRAPH_UNDIRECTED, modularity));
    }

    if (!membership) {
        igraph_vector_int_destroy(mymembership);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_le_community_to_membership
 * \brief Vertex membership from the leading eigenvector community structure.
 *
 * This function creates a membership vector from the
 * result of \ref igraph_community_leading_eigenvector().
 * It takes \c membership and performs \c steps merges,
 * according to the supplied \c merges matrix.
 *
 * \param merges The two-column matrix containing the merge operations.
 *    See \ref igraph_community_leading_eigenvector() for the
 *    detailed syntax. This is usually from the output of the
 *    leading eigenvector community structure detection routines.
 * \param steps The number of steps to make according to \c merges.
 * \param membership Initially the starting membership vector,
 *     on output the resulting membership vector, after performing \c steps merges.
 * \param csize Optionally the sizes of the communities are stored here,
 *     if this is not a null pointer, but an initialized vector.
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices.
 */
igraph_error_t igraph_le_community_to_membership(const igraph_matrix_int_t *merges,
                                                 igraph_integer_t steps,
                                                 igraph_vector_int_t *membership,
                                                 igraph_vector_int_t *csize) {

    igraph_integer_t no_of_nodes = igraph_vector_int_size(membership);
    igraph_vector_int_t fake_memb;
    igraph_integer_t components, i;

    if (no_of_nodes > 0) {
        components = igraph_vector_int_max(membership) + 1;
    } else {
        components = 0;
    }
    if (components > no_of_nodes) {
        IGRAPH_ERRORF("Invalid membership vector: number of components (%" IGRAPH_PRId ") must "
                      "not be greater than the number of nodes (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, components, no_of_nodes);
    }
    if (steps >= components) {
        IGRAPH_ERRORF("Number of steps (%" IGRAPH_PRId ") must be smaller than number of components (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, steps, components);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&fake_memb, components);

    /* Check membership vector */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*membership)[i] < 0) {
            IGRAPH_ERRORF("Invalid membership vector, negative ID found: %" IGRAPH_PRId ".", IGRAPH_EINVAL, VECTOR(*membership)[i]);
        }
        VECTOR(fake_memb)[ VECTOR(*membership)[i] ] += 1;
    }
    for (i = 0; i < components; i++) {
        if (VECTOR(fake_memb)[i] == 0) {
            /* Ideally the empty cluster's index would be reported.
               However, doing so would be confusing as some high-level interfaces
               use 1-based indexing, some 0-based. */
            IGRAPH_ERROR("Invalid membership vector, empty cluster found.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_community_to_membership(merges, components, steps, &fake_memb, 0));

    /* Ok, now we have the membership of the initial components,
       rewrite the original membership vector. */

    if (csize) {
        IGRAPH_CHECK(igraph_vector_int_resize(csize, components - steps));
        igraph_vector_int_null(csize);
    }

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*membership)[i] = VECTOR(fake_memb)[ VECTOR(*membership)[i] ];
        if (csize) {
            VECTOR(*csize)[ VECTOR(*membership)[i] ] += 1;
        }
    }

    igraph_vector_int_destroy(&fake_memb);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
