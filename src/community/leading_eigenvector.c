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
#include "igraph_memory.h"
#include "igraph_statusbar.h"
#include "igraph_structural.h"

#include "core/interruption.h"

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
    igraph_vector_t *idx;
    igraph_vector_t *idx2;
    igraph_adjlist_t *adjlist;
    igraph_inclist_t *inclist;
    igraph_vector_t *tmp;
    long int no_of_edges;
    igraph_vector_t *mymembership;
    long int comm;
    const igraph_vector_t *weights;
    const igraph_t *graph;
    igraph_vector_t *strength;
    igraph_real_t sumweights;
} igraph_i_community_leading_eigenvector_data_t;

static int igraph_i_community_leading_eigenvector(igraph_real_t *to,
                                                  const igraph_real_t *from,
                                                  int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_real_t ktx, ktx2;
    long int no_of_edges = data->no_of_edges;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        nlen = igraph_vector_int_size(neis);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int nei = (long int) VECTOR(*neis)[k];
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                to[j] += from[ (long int) VECTOR(*idx2)[nei] ];
                VECTOR(*tmp)[j] += 1;
            }
        }
    }

    /* Now calculate k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        long int degree = igraph_vector_int_size(neis);
        ktx += from[j] * degree;
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges / 2.0;
    ktx2 = ktx2 / no_of_edges / 2.0;

    /* Now calculate Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree = igraph_vector_int_size(neis);
        to[j] = to[j] - ktx * degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * degree;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static int igraph_i_community_leading_eigenvector2(igraph_real_t *to,
                                                   const igraph_real_t *from,
                                                   int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_adjlist_t *adjlist = data->adjlist;
    igraph_real_t ktx, ktx2;
    long int no_of_edges = data->no_of_edges;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        nlen = igraph_vector_int_size(neis);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int nei = (long int) VECTOR(*neis)[k];
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                long int fi = (long int) VECTOR(*idx2)[nei];
                if (fi < size) {
                    to[j] += from[fi];
                }
                VECTOR(*tmp)[j] += 1;
            }
        }
    }

    /* Now calculate k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size + 1; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        long int degree = igraph_vector_int_size(neis);
        if (j < size) {
            ktx += from[j] * degree;
        }
        ktx2 += degree;
    }
    ktx = ktx / no_of_edges / 2.0;
    ktx2 = ktx2 / no_of_edges / 2.0;

    /* Now calculate Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, oldid);
        igraph_real_t degree = igraph_vector_int_size(neis);
        to[j] = to[j] - ktx * degree;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * degree;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static int igraph_i_community_leading_eigenvector_weighted(igraph_real_t *to,
                                                           const igraph_real_t *from,
                                                           int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_inclist_t *inclist = data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *strength = data->strength;
    igraph_real_t sw = data->sumweights;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *inc = igraph_inclist_get(inclist, oldid);
        nlen = igraph_vector_int_size(inc);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int edge = (long int) VECTOR(*inc)[k];
            igraph_real_t w = VECTOR(*weights)[edge];
            long int nei = IGRAPH_OTHER(graph, edge, oldid);
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                to[j] += from[ (long int) VECTOR(*idx2)[nei] ] * w;
                VECTOR(*tmp)[j] += w;
            }
        }
    }

    /* k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        ktx += from[j] * str;
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;

    /* Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * str;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static int igraph_i_community_leading_eigenvector2_weighted(igraph_real_t *to,
                                                            const igraph_real_t *from,
                                                            int n, void *extra) {

    igraph_i_community_leading_eigenvector_data_t *data = extra;
    long int j, k, nlen, size = n;
    igraph_vector_t *idx = data->idx;
    igraph_vector_t *idx2 = data->idx2;
    igraph_vector_t *tmp = data->tmp;
    igraph_inclist_t *inclist = data->inclist;
    igraph_real_t ktx, ktx2;
    igraph_vector_t *mymembership = data->mymembership;
    long int comm = data->comm;
    const igraph_vector_t *weights = data->weights;
    const igraph_t *graph = data->graph;
    igraph_vector_t *strength = data->strength;
    igraph_real_t sw = data->sumweights;

    /* Ax */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_vector_int_t *inc = igraph_inclist_get(inclist, oldid);
        nlen = igraph_vector_int_size(inc);
        to[j] = 0.0;
        VECTOR(*tmp)[j] = 0.0;
        for (k = 0; k < nlen; k++) {
            long int edge = (long int) VECTOR(*inc)[k];
            igraph_real_t w = VECTOR(*weights)[edge];
            long int nei = IGRAPH_OTHER(graph, edge, oldid);
            long int neimemb = (long int) VECTOR(*mymembership)[nei];
            if (neimemb == comm) {
                long int fi = (long int) VECTOR(*idx2)[nei];
                if (fi < size) {
                    to[j] += from[fi] * w;
                }
                VECTOR(*tmp)[j] += w;
            }
        }
    }

    /* k^Tx/2m */
    ktx = 0.0; ktx2 = 0.0;
    for (j = 0; j < size + 1; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        if (j < size) {
            ktx += from[j] * str;
        }
        ktx2 += str;
    }
    ktx = ktx / sw / 2.0;
    ktx2 = ktx2 / sw / 2.0;

    /* Bx */
    for (j = 0; j < size; j++) {
        long int oldid = (long int) VECTOR(*idx)[j];
        igraph_real_t str = VECTOR(*strength)[oldid];
        to[j] = to[j] - ktx * str;
        VECTOR(*tmp)[j] = VECTOR(*tmp)[j] - ktx2 * str;
    }

    /* -d_ij summa l in G B_il */
    for (j = 0; j < size; j++) {
        to[j] -= VECTOR(*tmp)[j] * from[j];
    }

    return 0;
}

static void igraph_i_levc_free(igraph_vector_ptr_t *ptr) {
    long int i, n = igraph_vector_ptr_size(ptr);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*ptr)[i];
        if (v) {
            igraph_vector_destroy(v);
            igraph_free(v);
        }
    }
}

static void igraph_i_error_handler_none(const char *reason, const char *file,
                                        int line, int igraph_errno) {
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
 *
 * \param graph The undirected input graph.
 * \param weights The weights of the edges, or a null pointer for
 *    unweighted graphs.
 * \param merges The result of the algorithm, a matrix containing the
 *    information about the splits performed. The matrix is built in
 *    the opposite way however, it is like the result of an
 *    agglomerative algorithm. If at the end of the algorithm (after
 *    \p steps steps was done) there are <quote>p</quote> communities,
 *    then these are numbered from zero to <quote>p-1</quote>. The
 *    first line of the matrix contains the first <quote>merge</quote>
 *    (which is in reality the last split) of two communities into
 *    community <quote>p</quote>, the merge in the second line forms
 *    community <quote>p+1</quote>, etc. The matrix should be
 *    initialized before calling and will be resized as needed.
 *    This argument is ignored of it is \c NULL.
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
 * \param options The options for ARPACK. \c n is always
 *    overwritten. \c ncv is set to at least 4.
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
 *    that are calculated in each step of the algorithm, are stored here,
 *    in a pointer vector. Each eigenvector is stored in an
 *    \ref igraph_vector_t object. The user is responsible of
 *    deallocating the memory that belongs to the individual vectors,
 *    by calling first \ref igraph_vector_destroy(), and then
 *    \ref igraph_free() on them.
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
 *    calculation. If the callback returns a non-zero value, then the
 *    community finding algorithm stops. See the arguments passed to
 *    the callback at the documentation of \ref
 *    igraph_community_leading_eigenvector_callback_t.
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
int igraph_community_leading_eigenvector(const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_matrix_t *merges,
        igraph_vector_t *membership,
        igraph_integer_t steps,
        igraph_arpack_options_t *options,
        igraph_real_t *modularity,
        igraph_bool_t start,
        igraph_vector_t *eigenvalues,
        igraph_vector_ptr_t *eigenvectors,
        igraph_vector_t *history,
        igraph_community_leading_eigenvector_callback_t *callback,
        void *callback_extra) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_dqueue_t tosplit;
    igraph_vector_t idx, idx2, mymerges;
    igraph_vector_t strength, tmp;
    long int staken = 0;
    igraph_adjlist_t adjlist;
    igraph_inclist_t inclist;
    long int i, j, k, l;
    long int communities;
    igraph_vector_t vmembership, *mymembership = membership;
    igraph_i_community_leading_eigenvector_data_t extra;
    igraph_arpack_storage_t storage;
    igraph_real_t mod = 0;
    igraph_arpack_function_t *arpcb1 =
        weights ? igraph_i_community_leading_eigenvector_weighted :
        igraph_i_community_leading_eigenvector;
    igraph_arpack_function_t *arpcb2 =
        weights ? igraph_i_community_leading_eigenvector2_weighted :
        igraph_i_community_leading_eigenvector2;
    igraph_real_t sumweights = 0.0;

    if (weights && no_of_edges != igraph_vector_size(weights)) {
        IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }

    if (start && !membership) {
        IGRAPH_ERROR("Cannot start from given configuration if memberships "
                     "missing", IGRAPH_EINVAL);
    }

    if (start && membership &&
        igraph_vector_size(membership) != no_of_nodes) {
        IGRAPH_ERROR("Wrong length for vector of predefined memberships",
                     IGRAPH_EINVAL);
    }

    if (start && membership && igraph_vector_max(membership) >= no_of_nodes) {
        IGRAPH_WARNING("Too many communities in membership start vector");
    }

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("This method was developed for undirected graphs");
    }

    if (steps < 0 || steps > no_of_nodes - 1) {
        steps = (igraph_integer_t) no_of_nodes - 1;
    }

    if (!membership) {
        mymembership = &vmembership;
        IGRAPH_VECTOR_INIT_FINALLY(mymembership, 0);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&mymerges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&mymerges, steps * 2));
    IGRAPH_VECTOR_INIT_FINALLY(&idx, 0);
    if (eigenvalues)  {
        igraph_vector_clear(eigenvalues);
    }
    if (eigenvectors) {
        igraph_vector_ptr_clear(eigenvectors);
        IGRAPH_FINALLY(igraph_i_levc_free, eigenvectors);
    }

    IGRAPH_STATUS("Starting leading eigenvector method.\n", 0);

    if (!start) {
        /* Calculate the weakly connected components in the graph and use them as
         * an initial split */
        IGRAPH_CHECK(igraph_clusters(graph, mymembership, &idx, 0, IGRAPH_WEAK));
        communities = igraph_vector_size(&idx);
        IGRAPH_STATUSF(("Starting from %li component(s).\n", 0, communities));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history,
                                                 IGRAPH_LEVC_HIST_START_FULL));
        }
    } else {
        /* Just create the idx vector for the given membership vector */
        communities = (long int) igraph_vector_max(mymembership) + 1;
        IGRAPH_STATUSF(("Starting from given membership vector with %li "
                        "communities.\n", 0, communities));
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history,
                                                 IGRAPH_LEVC_HIST_START_GIVEN));
            IGRAPH_CHECK(igraph_vector_push_back(history, communities));
        }
        IGRAPH_CHECK(igraph_vector_resize(&idx, communities));
        igraph_vector_null(&idx);
        for (i = 0; i < no_of_nodes; i++) {
            int t = (int) VECTOR(*mymembership)[i];
            VECTOR(idx)[t] += 1;
        }
    }

    IGRAPH_DQUEUE_INIT_FINALLY(&tosplit, 100);
    for (i = 0; i < communities; i++) {
        if (VECTOR(idx)[i] > 2) {
            igraph_dqueue_push(&tosplit, i);
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
            igraph_vector_t *v = IGRAPH_CALLOC(1, igraph_vector_t);
            if (!v) {
                IGRAPH_ERROR("Cannot do leading eigenvector community detection",
                             IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, v);
            IGRAPH_VECTOR_INIT_FINALLY(v, 0);
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eigenvectors, v));
            IGRAPH_FINALLY_CLEAN(2);
        }
        if (history) {
            IGRAPH_CHECK(igraph_vector_push_back(history, IGRAPH_LEVC_HIST_SPLIT));
            IGRAPH_CHECK(igraph_vector_push_back(history, i - 1));
        }
    }
    staken = communities - 1;

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_resize(&idx, no_of_nodes));
    igraph_vector_null(&idx);
    IGRAPH_VECTOR_INIT_FINALLY(&idx2, no_of_nodes);
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

    options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
    options->start = 0;
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

    while (!igraph_dqueue_empty(&tosplit) && staken < steps) {
        long int comm = (long int) igraph_dqueue_pop_back(&tosplit);
        /* depth first search */
        long int size = 0;
        igraph_real_t tmpev;

        IGRAPH_STATUSF(("Trying to split community %li... ", 0, comm));
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

        /* We solve two eigenproblems, one for the original modularity
           matrix, and one for the modularity matrix after deleting the
           last row and last column from it. This is a trick to find
           multiple leading eigenvalues, because ARPACK is sometimes
           unstable when the first two eigenvalues are requested, but it
           does much better for the single principal eigenvalue. */

        /* We start with the smaller eigenproblem. */

        options->n = (int) size - 1;
        options->info = 0;
        options->nev = 1;
        options->ldv = 0;
        options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */
        options->nconv = 0;
        options->lworkl = 0;        /* we surely have enough space */
        extra.comm = comm;

        /* We try calling the solver twice, once from a random starting
           point, once from a fixed one. This is because for some hard
           cases it tends to fail. We need to suppress error handling for
           the first call. */
        {
            int i;
            igraph_error_handler_t *errh =
                igraph_set_error_handler(igraph_i_error_handler_none);
            igraph_warning_handler_t *warnh =
                igraph_set_warning_handler(igraph_warning_handler_ignore);
            igraph_arpack_rssolve(arpcb2, &extra, options, &storage,
                                  /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            igraph_set_warning_handler(warnh);
            if (options->nconv < 1) {
                /* Call again from a fixed starting point. Note that we cannot use a
                 * fixed all-1 starting vector as sometimes ARPACK would return a
                 * 'starting vector is zero' error -- this is of course not true but
                 * it's a result of ARPACK >= 3.6.3 trying to force the starting vector
                 * into the range of OP (i.e. the matrix being solved). The initial
                 * vector we use here seems to work, but I have no theoretical argument
                 * for its usage; it just happens to work. */
                options->start = 1;
                options->info = 0;
                options->ncv = 0;
                options->lworkl = 0;    /* we surely have enough space */
                for (i = 0; i < options->n ; i++) {
                    storage.resid[i] = i % 2 ? 1 : -1;
                }
                IGRAPH_CHECK(igraph_arpack_rssolve(arpcb2, &extra, options, &storage,
                                                   /*values=*/ 0, /*vectors=*/ 0));
                options->start = 0;
            }
        }

        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }

        tmpev = storage.d[0];

        /* Now we do the original eigenproblem, again, twice if needed */

        options->n = (int) size;
        options->info = 0;
        options->nev = 1;
        options->ldv = 0;
        options->nconv = 0;
        options->lworkl = 0;    /* we surely have enough space */
        options->ncv = 0;   /* 0 means "automatic" in igraph_arpack_rssolve */

        {
            int i;
            igraph_error_handler_t *errh =
                igraph_set_error_handler(igraph_i_error_handler_none);
            igraph_arpack_rssolve(arpcb1, &extra, options, &storage,
                                  /*values=*/ 0, /*vectors=*/ 0);
            igraph_set_error_handler(errh);
            if (options->nconv < 1) {
                /* Call again from a fixed starting point. See the comment a few lines
                 * above about the exact choice of this starting vector */
                options->start = 1;
                options->info = 0;
                options->ncv = 0;
                options->lworkl = 0;    /* we surely have enough space */
                for (i = 0; i < options->n; i++) {
                    storage.resid[i] = i % 2 ? 1 : -1;
                }
                IGRAPH_CHECK(igraph_arpack_rssolve(arpcb1, &extra, options, &storage,
                                                   /*values=*/ 0, /*vectors=*/ 0));
                options->start = 0;
            }
        }

        if (options->nconv < 1) {
            IGRAPH_ERROR("ARPACK did not converge", IGRAPH_ARPACK_FAILED);
        }

        /* Ok, we have the leading eigenvector of the modularity matrix*/

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
            int ret;
            igraph_vector_view(&vv, storage.v, size);
            ret = callback(mymembership, comm, storage.d[0], &vv,
                           arpcb1, &extra, callback_extra);
            if (ret) {
                break;
            }
        }

        if (eigenvalues) {
            IGRAPH_CHECK(igraph_vector_push_back(eigenvalues, storage.d[0]));
        }

        if (eigenvectors) {
            igraph_vector_t *v = IGRAPH_CALLOC(1, igraph_vector_t);
            if (!v) {
                IGRAPH_ERROR("Cannot do leading eigenvector community detection",
                             IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, v);
            IGRAPH_VECTOR_INIT_FINALLY(v, size);
            for (i = 0; i < size; i++) {
                VECTOR(*v)[i] = storage.v[i];
            }
            IGRAPH_CHECK(igraph_vector_ptr_push_back(eigenvectors, v));
            IGRAPH_FINALLY_CLEAN(2);
        }

        if (storage.d[0] <= 0) {
            IGRAPH_STATUS("no split.\n", 0);
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        /* Check for multiple leading eigenvalues */

        if (fabs(storage.d[0] - tmpev) < 1e-8) {
            IGRAPH_STATUS("multiple principal eigenvalue, no split.\n", 0);
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
            IGRAPH_STATUS("no split.\n", 0);
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
            IGRAPH_STATUS("no modularity increase, no split.\n", 0);
            if (history) {
                IGRAPH_CHECK(igraph_vector_push_back(history,
                                                     IGRAPH_LEVC_HIST_FAILED));
                IGRAPH_CHECK(igraph_vector_push_back(history, comm));
            }
            continue;
        }

        communities++;
        IGRAPH_STATUS("split.\n", 0);

        /* Rewrite the mymembership vector */
        for (j = 0; j < size; j++) {
            if (storage.v[j] < 0) {
                long int oldid = (long int) VECTOR(idx)[j];
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
            IGRAPH_CHECK(igraph_dqueue_push(&tosplit, communities - 1));
        }
        if (size - l > 1) {
            IGRAPH_CHECK(igraph_dqueue_push(&tosplit, comm));
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
    igraph_dqueue_destroy(&tosplit);
    igraph_vector_destroy(&tmp);
    igraph_vector_destroy(&idx2);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_STATUS("Done.\n", 0);

    /* reform the mymerges vector */
    if (merges) {
        igraph_vector_null(&idx);
        l = igraph_vector_size(&mymerges);
        k = communities;
        j = 0;
        IGRAPH_CHECK(igraph_matrix_resize(merges, l / 2, 2));
        for (i = l; i > 0; i -= 2) {
            long int from = (long int) VECTOR(mymerges)[i - 1];
            long int to = (long int) VECTOR(mymerges)[i - 2];
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

    if (eigenvectors) {
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&idx);
    igraph_vector_destroy(&mymerges);
    IGRAPH_FINALLY_CLEAN(2);

    if (modularity) {
      IGRAPH_CHECK(igraph_modularity(graph, mymembership, weights,
                                    /* resolution */ 1,
                                    /* only undirected */ 0, modularity));
    }

    if (!membership) {
        igraph_vector_destroy(mymembership);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_le_community_to_membership
 * Vertex membership from the leading eigenvector community structure
 *
 * This function creates a membership vector from the
 * result of \ref igraph_community_leading_eigenvector(),
 * It takes \c membership
 * and performs \c steps merges, according to the supplied
 * \c merges matrix.
 * \param merges The two-column matrix containing the merge
 *    operations. See \ref igraph_community_walktrap() for the
 *    detailed syntax. This is usually from the output of the
 *    leading eigenvector community structure detection routines.
 * \param steps The number of steps to make according to \c merges.
 * \param membership Initially the starting membership vector,
 *     on output the resulting membership vector, after performing \c steps merges.
 * \param csize Optionally the sizes of the communities is stored here,
 *     if this is not a null pointer, but an initialized vector.
 * \return Error code.
 *
 * Time complexity: O(|V|), the number of vertices.
 */
int igraph_le_community_to_membership(const igraph_matrix_t *merges,
                                      igraph_integer_t steps,
                                      igraph_vector_t *membership,
                                      igraph_vector_t *csize) {

    long int no_of_nodes = igraph_vector_size(membership);
    igraph_vector_t fake_memb;
    long int components, i;

    if (no_of_nodes > 0) {
        components = (long int) igraph_vector_max(membership) + 1;
    } else {
        components = 0;
    }
    if (components > no_of_nodes) {
        IGRAPH_ERRORF("Invalid membership vector: number of components (%ld) must "
         "not be greater than the number of nodes (%ld).", IGRAPH_EINVAL, components, no_of_nodes);
    }
    if (steps >= components) {
        IGRAPH_ERRORF("Number of steps (%" IGRAPH_PRId ") must be smaller than number of components (%ld).",
                      IGRAPH_EINVAL, steps, components);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&fake_memb, components);

    /* Check membership vector */
    for (i = 0; i < no_of_nodes; i++) {
        if (VECTOR(*membership)[i] < 0) {
            IGRAPH_ERRORF("Invalid membership vector, negative ID found: %g.", IGRAPH_EINVAL, VECTOR(*membership)[i]);
        }
        VECTOR(fake_memb)[ (long int) VECTOR(*membership)[i] ] += 1;
    }
    for (i = 0; i < components; i++) {
        if (VECTOR(fake_memb)[i] == 0) {
            /* Ideally the empty cluster's index would be reported.
               However, doing so would be confusing as some high-level interfaces
               use 1-based indexing, some 0-based. */
            IGRAPH_ERROR("Invalid membership vector, empty cluster found.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_community_to_membership(merges, (igraph_integer_t)
                 components, steps,
                 &fake_memb, 0));

    /* Ok, now we have the membership of the initial components,
       rewrite the original membership vector. */

    if (csize) {
        IGRAPH_CHECK(igraph_vector_resize(csize, components - steps));
        igraph_vector_null(csize);
    }

    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(*membership)[i] = VECTOR(fake_memb)[ (long int) VECTOR(*membership)[i] ];
        if (csize) {
            VECTOR(*csize)[ (long int) VECTOR(*membership)[i] ] += 1;
        }
    }

    igraph_vector_destroy(&fake_memb);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
