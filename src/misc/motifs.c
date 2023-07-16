/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_motifs.h"

#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_nongraph.h"
#include "igraph_stack.h"

#include "core/interruption.h"
#include "isomorphism/isoclasses.h"
#include "graph/internal.h"

/**
 * Callback function for igraph_motifs_randesu that counts the motifs by
 * isomorphism class in a histogram.
 */
static igraph_error_t igraph_i_motifs_randesu_update_hist(
        const igraph_t *graph,
        igraph_vector_int_t *vids, igraph_integer_t isoclass, void* extra) {
    igraph_vector_t *hist = (igraph_vector_t*)extra;
    IGRAPH_UNUSED(graph); IGRAPH_UNUSED(vids);
    VECTOR(*hist)[isoclass]++;
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_motifs_randesu
 * \brief Count the number of motifs in a graph.
 *
 * </para><para>
 * Motifs are small weakly connected induced subgraphs of a given structure in a
 * graph. It is argued that the motif profile (i.e. the number of
 * different motifs in the graph) is characteristic for different
 * types of networks and network function is related to the motifs in
 * the graph.
 *
 * </para><para>
 * This function is able to find directed motifs of sizes three
 * and four and undirected motifs of sizes three to six
 * (i.e. the number of different subgraphs with three to six
 * vertices in the network).
 *
 * </para><para>
 * In a big network the total number of motifs can be very large, so
 * it takes a lot of time to find all of them. In this case, a sampling
 * method can be used. This function is capable of doing sampling via the
 * \p cut_prob argument. This argument gives the probability that
 * a branch of the motif search tree will not be explored. See
 * S. Wernicke and F. Rasche: FANMOD: a tool for fast network motif
 * detection, Bioinformatics 22(9), 1152--1153, 2006 for details.
 * https://doi.org/10.1093/bioinformatics/btl038
 *
 * </para><para>
 * Set the \p cut_prob argument to a zero vector for finding all
 * motifs.
 *
 * </para><para>
 * Directed motifs will be counted in directed graphs and undirected
 * motifs in undirected graphs.
 *
 * \param graph The graph to find the motifs in.
 * \param hist The result of the computation, it gives the number of
 *        motifs found for each isomorphism class. See
 *        \ref igraph_isoclass() for help about isomorphism classes.
 *        Note that this function does \em not count isomorphism
 *        classes that are not connected and will report NaN (more
 *        precisely \c IGRAPH_NAN) for them.
 * \param size The size of the motifs to search for. For directed graphs,
 *        only 3 and 4 are implemented, for undirected, 3 to 6.
 *        The limitation is not in the motif finding code, but the graph
 *        isomorphism code.
 * \param cut_prob Vector of probabilities for cutting the search tree
 *        at a given level. The first element is the first level, etc.
 *        Supply all zeros here (of length \p size) to find all motifs
 *        in a graph.
 * \return Error code.
 *
 * \sa \ref igraph_motifs_randesu_estimate() for estimating the number
 * of motifs in a graph, this can help to set the \p cut_prob
 * parameter; \ref igraph_motifs_randesu_no() to calculate the total
 * number of motifs of a given size in a graph;
 * \ref igraph_motifs_randesu_callback() for calling a callback function
 * for every motif found; \ref igraph_subisomorphic_lad() for finding
 * subgraphs on more than 4 (directed) or 6 (undirected) vertices;
 * \ref igraph_graph_count() to find the number of graph on a given
 * number of vertices, i.e. the length of the \p hist vector.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/igraph_motifs_randesu.c
 */
igraph_error_t igraph_motifs_randesu(const igraph_t *graph, igraph_vector_t *hist,
                          igraph_integer_t size, const igraph_vector_t *cut_prob) {
    igraph_bool_t directed = igraph_is_directed(graph);
    igraph_integer_t histlen;

    if (directed) {
        switch (size) {
        case 3:
            histlen = 16;
            break;
        case 4:
            histlen = 218;
            break;
        default:
            IGRAPH_ERROR("In directed graphs, only 3 and 4 vertex motifs are supported.",
                         IGRAPH_UNIMPLEMENTED);
        }
    } else {
        switch (size) {
        case 3:
            histlen = 4;
            break;
        case 4:
            histlen = 11;
            break;
        case 5:
            histlen = 34;
            break;
        case 6:
            histlen = 156;
            break;
        default:
            IGRAPH_ERROR("In undirected graphs, only 3 to 6 vertex motifs are supported.",
                         IGRAPH_UNIMPLEMENTED);
        }
    }

    if (igraph_vector_size(cut_prob) != size) {
        IGRAPH_ERRORF("Cut probability vector size (%" IGRAPH_PRId ") must agree with motif size (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(cut_prob), size);
    }

    IGRAPH_CHECK(igraph_vector_resize(hist, histlen));
    igraph_vector_null(hist);

    IGRAPH_CHECK(igraph_motifs_randesu_callback(graph, size, cut_prob,
                 &igraph_i_motifs_randesu_update_hist, hist));

    if (size == 3) {
        if (directed) {
            VECTOR(*hist)[0] = VECTOR(*hist)[1] = VECTOR(*hist)[3] = IGRAPH_NAN;
        } else {
            VECTOR(*hist)[0] = VECTOR(*hist)[1] = IGRAPH_NAN;
        }
    } else if (size == 4) {
        if (directed) {
            const int not_connected[] = { 0, 1, 2, 4, 5, 6, 9, 10, 11, 15, 22, 23, 27,
                                          28, 33, 34, 39, 62, 120 };
            size_t i, n = sizeof(not_connected) / sizeof(not_connected[0]);
            for (i = 0; i < n; i++) {
                VECTOR(*hist)[not_connected[i]] = IGRAPH_NAN;
            }
        } else {
            VECTOR(*hist)[0] = VECTOR(*hist)[1] = VECTOR(*hist)[2] =
                    VECTOR(*hist)[3] = VECTOR(*hist)[5] = IGRAPH_NAN;
        }
    } else if (size == 5) {
        /* undirected only */
        const int not_connected[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 19 };
        size_t i, n = sizeof(not_connected) / sizeof(int);
        for (i = 0; i < n; i++) {
            VECTOR(*hist)[not_connected[i]] = IGRAPH_NAN;
        }
    } else if (size == 6) {
        /* undirected only */
        const int not_connected[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                     16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                     30, 31, 32, 33, 35, 38, 44, 50, 51, 54, 74, 77, 89, 120};
        size_t i, n = sizeof(not_connected) / sizeof(int);
        for (i = 0; i < n; i++) {
            VECTOR(*hist)[not_connected[i]] = IGRAPH_NAN;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_motifs_randesu_callback
 * \brief Finds motifs in a graph and calls a function for each of them.
 *
 * </para><para>
 * Similarly to \ref igraph_motifs_randesu(), this function is able to find
 * directed motifs of sizes three and four and undirected motifs of sizes
 * three to six (i.e. the number of different subgraphs with three to six
 * vertices in the network). However, instead of
 * counting them, the function will call a callback function for each motif
 * found to allow further tests or post-processing.
 *
 * </para><para>
 * The \p cut_prob argument also allows sampling the motifs, just like for
 * \ref igraph_motifs_randesu(). Set the \p cut_prob argument to a zero vector
 * for finding all motifs.
 *
 * \param graph The graph to find the motifs in.
 * \param size The size of the motifs to search for. Only three and
 *        four are implemented currently. The limitation is not in the
 *        motif finding code, but the graph isomorphism code.
 * \param cut_prob Vector of probabilities for cutting the search tree
 *        at a given level. The first element is the first level, etc.
 *        Supply all zeros here (of length \c size) to find all motifs
 *        in a graph.
 * \param callback A pointer to a function of type \ref igraph_motifs_handler_t.
 *        This function will be called whenever a new motif is found.
 * \param extra Extra argument to pass to the callback function.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/igraph_motifs_randesu.c
 */

igraph_error_t igraph_motifs_randesu_callback(const igraph_t *graph, igraph_integer_t size,
                                   const igraph_vector_t *cut_prob, igraph_motifs_handler_t *callback,
                                   void* extra) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_adjlist_t allneis, alloutneis;
    igraph_vector_int_t *neis;
    igraph_integer_t father;
    igraph_integer_t i, j, s;
    igraph_integer_t motifs = 0;
    IGRAPH_UNUSED(motifs);    /* We mark it as unused to prevent warnings about unused-but-set-variables. */

    igraph_vector_int_t vids;     /* this is G */
    igraph_vector_int_t adjverts; /* this is V_E */
    igraph_stack_int_t stack;     /* this is S */
    igraph_integer_t *added;
    char *subg;

    const unsigned int *arr_idx, *arr_code;
    unsigned int code = 0;
    unsigned int mul, idx;

    igraph_bool_t terminate = false;

    if (igraph_is_directed(graph)) {
        switch (size) {
        case 3:
            arr_idx = igraph_i_isoclass_3_idx;
            arr_code = igraph_i_isoclass2_3;
            mul = 3;
            break;
        case 4:
            arr_idx = igraph_i_isoclass_4_idx;
            arr_code = igraph_i_isoclass2_4;
            mul = 4;
            break;
        default:
            IGRAPH_ERROR("In directed graphs, only 3 and 4 vertex motifs are supported.",
                         IGRAPH_UNIMPLEMENTED);
        }
    } else {
        switch (size) {
        case 3:
            arr_idx = igraph_i_isoclass_3u_idx;
            arr_code = igraph_i_isoclass2_3u;
            mul = 3;
            break;
        case 4:
            arr_idx = igraph_i_isoclass_4u_idx;
            arr_code = igraph_i_isoclass2_4u;
            mul = 4;
            break;
        case 5:
            arr_idx = igraph_i_isoclass_5u_idx;
            arr_code = igraph_i_isoclass2_5u;
            mul = 5;
            break;
        case 6:
            arr_idx = igraph_i_isoclass_6u_idx;
            arr_code = igraph_i_isoclass2_6u;
            mul = 6;
            break;
        default:
            IGRAPH_ERROR("In undirected graphs, only 3 to 6 vertex motifs are supported.",
                         IGRAPH_UNIMPLEMENTED);
        }
    }

    if (igraph_vector_size(cut_prob) != size) {
        IGRAPH_ERRORF("Cut probability vector size (%" IGRAPH_PRId ") must agree with motif size (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(cut_prob), size);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Insufficient memory to find motifs.");
    IGRAPH_FINALLY(igraph_free, added);

    subg = IGRAPH_CALLOC(no_of_nodes, char);
    IGRAPH_CHECK_OOM(subg, "Insufficient memory to find motifs.");
    IGRAPH_FINALLY(igraph_free, subg);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &alloutneis, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &alloutneis);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&adjverts, 0);
    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);

    RNG_BEGIN();

    for (father = 0; father < no_of_nodes; father++) {
        igraph_integer_t level;

        IGRAPH_ALLOW_INTERRUPTION();

        if (VECTOR(*cut_prob)[0] == 1 || RNG_UNIF01() < VECTOR(*cut_prob)[0]) {
            continue;
        }

        /* init G */
        igraph_vector_int_clear(&vids); level = 0;
        IGRAPH_CHECK(igraph_vector_int_push_back(&vids, father));
        subg[father] = 1; added[father] += 1; level += 1;

        /* init V_E */
        igraph_vector_int_clear(&adjverts);
        neis = igraph_adjlist_get(&allneis, father);
        s = igraph_vector_int_size(neis);
        for (i = 0; i < s; i++) {
            igraph_integer_t nei = VECTOR(*neis)[i];
            if (!added[nei] && nei > father) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei));
                IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, father));
            }
            added[nei] += 1;
        }

        /* init S */
        igraph_stack_int_clear(&stack);

        while (level > 1 || !igraph_vector_int_empty(&adjverts)) {
            igraph_real_t cp = VECTOR(*cut_prob)[level];

            if (level == size - 1) {
                s = igraph_vector_int_size(&adjverts) / 2;
                for (i = 0; i < s; i++) {
                    igraph_integer_t k, s2;
                    igraph_integer_t last;
                    igraph_error_t ret;

                    if (cp != 0 && RNG_UNIF01() < cp) {
                        continue;
                    }
                    motifs += 1;

                    last = VECTOR(adjverts)[2 * i];
                    IGRAPH_CHECK(igraph_vector_int_push_back(&vids, last));
                    subg[last] = (char) size;

                    code = 0; idx = 0;
                    for (k = 0; k < size; k++) {
                        igraph_integer_t from = VECTOR(vids)[k];
                        neis = igraph_adjlist_get(&alloutneis, from);
                        s2 = igraph_vector_int_size(neis);
                        for (j = 0; j < s2; j++) {
                            igraph_integer_t nei = VECTOR(*neis)[j];
                            if (subg[nei] && k != subg[nei] - 1) {
                                idx = (unsigned char) (mul * k + (subg[nei] - 1));
                                code |= arr_idx[idx];
                            }
                        }
                    }

                    IGRAPH_CHECK_CALLBACK(
                        callback(graph, &vids, arr_code[code], extra),
                        &ret
                    );

                    if (ret == IGRAPH_STOP) {
                        terminate = true;
                        break;
                    }

                    igraph_vector_int_pop_back(&vids);
                    subg[last] = 0;
                }
            }

            /* did the callback function asked us to terminate the search? */
            if (terminate) {
                break;
            }

            /* can we step down? */
            if (level < size - 1 &&
                !igraph_vector_int_empty(&adjverts)) {
                /* we might step down */
                igraph_integer_t neifather = igraph_vector_int_pop_back(&adjverts);
                igraph_integer_t nei = igraph_vector_int_pop_back(&adjverts);

                if (cp == 0 || RNG_UNIF01() > cp) {
                    /* yes, step down */
                    IGRAPH_CHECK(igraph_vector_int_push_back(&vids, nei));
                    subg[nei] = (char) level + 1; added[nei] += 1; level += 1;

                    IGRAPH_CHECK(igraph_stack_int_push(&stack, neifather));
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, nei));
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, level));

                    neis = igraph_adjlist_get(&allneis, nei);
                    s = igraph_vector_int_size(neis);
                    for (i = 0; i < s; i++) {
                        igraph_integer_t nei2 = VECTOR(*neis)[i];
                        if (!added[nei2] && nei2 > father) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei2));
                            IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei));
                        }
                        added[nei2] += 1;
                    }
                }
            } else {
                /* no, step back */
                igraph_integer_t nei, neifather;
                while (!igraph_stack_int_empty(&stack) &&
                       level == igraph_stack_int_top(&stack) - 1) {
                    igraph_stack_int_pop(&stack);
                    nei = igraph_stack_int_pop(&stack);
                    neifather = igraph_stack_int_pop(&stack);
                    igraph_vector_int_push_back(&adjverts, nei);
                    igraph_vector_int_push_back(&adjverts, neifather);
                }

                nei = igraph_vector_int_pop_back(&vids);
                subg[nei] = 0; added[nei] -= 1; level -= 1;
                neis = igraph_adjlist_get(&allneis, nei);
                s = igraph_vector_int_size(neis);
                for (i = 0; i < s; i++) {
                    added[ VECTOR(*neis)[i] ] -= 1;
                }
                while (!igraph_vector_int_empty(&adjverts) &&
                       igraph_vector_int_tail(&adjverts) == nei) {
                    igraph_vector_int_pop_back(&adjverts);
                    igraph_vector_int_pop_back(&adjverts);
                }
            }

        } /* while */

        /* did the callback function asked us to terminate the search? */
        if (terminate) {
            break;
        }

        /* clear the added vector */
        added[father] -= 1;
        subg[father] = 0;
        neis = igraph_adjlist_get(&allneis, father);
        s = igraph_vector_int_size(neis);
        for (i = 0; i < s; i++) {
            added[ VECTOR(*neis)[i] ] -= 1;
        }

    } /* for father */

    RNG_END();

    IGRAPH_FREE(added);
    IGRAPH_FREE(subg);
    igraph_vector_int_destroy(&vids);
    igraph_vector_int_destroy(&adjverts);
    igraph_adjlist_destroy(&alloutneis);
    igraph_adjlist_destroy(&allneis);
    igraph_stack_int_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_motifs_randesu_estimate
 * \brief Estimate the total number of motifs in a graph.
 *
 * This function estimates the total number of (weakly) connected induced
 * subgraphs on \p size vertices. For example, an undirected complete graph
 * on \c n vertices will have one motif of size \c n, and \c n motifs
 * of \p size <code>n - 1</code>. As another example, one triangle
 * and a separate vertex will have zero motifs of size four.
 *
 * </para><para>
 * This function is useful for large graphs for which it is not
 * feasible to count all connected subgraphs, as there are too
 * many of them.
 *
 * </para><para>
 * The estimate is made by taking a sample of vertices and counting all
 * connected subgraphs in which these vertices are included. There is also
 * a \p cut_prob parameter which gives the probabilities to cut a branch of
 * the search tree.
 *
 * \param graph The graph object to study.
 * \param est Pointer to an integer, the result will be stored here.
 * \param size The size of the subgraphs to look for.
 * \param cut_prob Vector giving the probabilities to cut a branch of
 *        the search tree and omit counting the motifs in that branch.
 *        It contains a probability for each level. Supply \p size
 *        zeros here to count all the motifs in the sample.
 * \param sample_size The number of vertices to use as the
 *        sample. This parameter is only used if the \p parsample
 *        argument is a null pointer.
 * \param parsample Either pointer to an initialized vector or a null
 *        pointer. If a vector then the vertex IDs in the vector are
 *        used as a sample. If a null pointer then the \p sample_size
 *        argument is used to create a sample of vertices drawn with
 *        uniform probability.
 * \return Error code.
 *
 * \sa \ref igraph_motifs_randesu(), \ref igraph_motifs_randesu_no().
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_motifs_randesu_estimate(const igraph_t *graph, igraph_integer_t *est,
                                   igraph_integer_t size, const igraph_vector_t *cut_prob,
                                   igraph_integer_t sample_size,
                                   const igraph_vector_int_t *parsample) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t neis;

    igraph_vector_int_t vids;     /* this is G */
    igraph_vector_int_t adjverts; /* this is V_E */
    igraph_stack_int_t stack;     /* this is S */
    igraph_integer_t *added;
    igraph_vector_int_t *sample;
    igraph_integer_t sam;
    igraph_integer_t i;

    if (size < 3) {
        IGRAPH_ERRORF("Motif size must be at least 3, received %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, size);
    }

    if (igraph_vector_size(cut_prob) != size) {
        IGRAPH_ERRORF("Cut probability vector size (%" IGRAPH_PRId ") must agree with motif size (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(cut_prob), size);
    }

    if (parsample && !igraph_vector_int_isininterval(parsample, 0, no_of_nodes-1)) {
        IGRAPH_ERROR("Sample vertex ID out of range.", IGRAPH_EINVVID);
    }

    if (no_of_nodes == 0) {
        *est = 0;
        return IGRAPH_SUCCESS;
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Insufficient memory to count motifs.");
    IGRAPH_FINALLY(igraph_free, added);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&adjverts, 0);
    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    if (parsample == NULL) {
        sample = IGRAPH_CALLOC(1, igraph_vector_int_t);
        IGRAPH_CHECK_OOM(sample, "Insufficient memory to count motifs.");
        IGRAPH_FINALLY(igraph_free, sample);
        IGRAPH_VECTOR_INT_INIT_FINALLY(sample, 0);
        IGRAPH_CHECK(igraph_random_sample(sample, 0, no_of_nodes - 1, sample_size));
    } else {
        sample = (igraph_vector_int_t*) parsample;
        sample_size = igraph_vector_int_size(sample);
    }

    *est = 0;

    RNG_BEGIN();

    for (sam = 0; sam < sample_size; sam++) {
        igraph_integer_t father = VECTOR(*sample)[sam];
        igraph_integer_t level, s;

        IGRAPH_ALLOW_INTERRUPTION();

        if (VECTOR(*cut_prob)[0] == 1 ||
            RNG_UNIF01() < VECTOR(*cut_prob)[0]) {
            continue;
        }

        /* init G */
        igraph_vector_int_clear(&vids); level = 0;
        IGRAPH_CHECK(igraph_vector_int_push_back(&vids, father));
        added[father] += 1; level += 1;

        /* init V_E */
        igraph_vector_int_clear(&adjverts);
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
        s = igraph_vector_int_size(&neis);
        for (i = 0; i < s; i++) {
            igraph_integer_t nei = VECTOR(neis)[i];
            if (!added[nei] && nei > father) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei));
                IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, father));
            }
            added[nei] += 1;
        }

        /* init S */
        igraph_stack_int_clear(&stack);

        while (level > 1 || !igraph_vector_int_empty(&adjverts)) {
            igraph_real_t cp = VECTOR(*cut_prob)[level];

            if (level == size - 1) {
                s = igraph_vector_int_size(&adjverts) / 2;
                for (i = 0; i < s; i++) {
                    if (cp != 0 && RNG_UNIF01() < cp) {
                        continue;
                    }
                    (*est) += 1;
                }
            }

            if (level < size - 1 &&
                !igraph_vector_int_empty(&adjverts)) {
                /* We might step down */
                igraph_integer_t neifather = igraph_vector_int_pop_back(&adjverts);
                igraph_integer_t nei = igraph_vector_int_pop_back(&adjverts);

                if (cp == 0 || RNG_UNIF01() > cp) {
                    /* Yes, step down */
                    IGRAPH_CHECK(igraph_vector_int_push_back(&vids, nei));
                    added[nei] += 1; level += 1;

                    IGRAPH_CHECK(igraph_stack_int_push(&stack, neifather));
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, nei));
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, level));

                    IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
                    s = igraph_vector_int_size(&neis);
                    for (i = 0; i < s; i++) {
                        igraph_integer_t nei2 = VECTOR(neis)[i];
                        if (!added[nei2] && nei2 > father) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei2));
                            IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei));
                        }
                        added[nei2] += 1;
                    }
                }
            } else {
                /* no, step back */
                igraph_integer_t nei, neifather;
                while (!igraph_stack_int_empty(&stack) &&
                       level == igraph_stack_int_top(&stack) - 1) {
                    igraph_stack_int_pop(&stack);
                    nei = igraph_stack_int_pop(&stack);
                    neifather = igraph_stack_int_pop(&stack);
                    igraph_vector_int_push_back(&adjverts, nei);
                    igraph_vector_int_push_back(&adjverts, neifather);
                }

                nei = igraph_vector_int_pop_back(&vids);
                added[nei] -= 1; level -= 1;
                IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
                s = igraph_vector_int_size(&neis);
                for (i = 0; i < s; i++) {
                    added[ VECTOR(neis)[i] ] -= 1;
                }
                while (!igraph_vector_int_empty(&adjverts) &&
                       igraph_vector_int_tail(&adjverts) == nei) {
                    igraph_vector_int_pop_back(&adjverts);
                    igraph_vector_int_pop_back(&adjverts);
                }
            }

        } /* while */

        /* clear the added vector */
        added[father] -= 1;
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
        s = igraph_vector_int_size(&neis);
        for (i = 0; i < s; i++) {
            added[ VECTOR(neis)[i] ] -= 1;
        }

    } /* for father */

    RNG_END();

    (*est) *= ((igraph_real_t) no_of_nodes / sample_size);

    if (parsample == 0) {
        igraph_vector_int_destroy(sample);
        IGRAPH_FREE(sample);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FREE(added);
    igraph_vector_int_destroy(&vids);
    igraph_vector_int_destroy(&adjverts);
    igraph_stack_int_destroy(&stack);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(5);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_motifs_randesu_no
 * \brief Count the total number of motifs in a graph.
 *
 * This function counts the total number of (weakly) connected
 * induced subgraphs on \p size vertices, without assigning isomorphism
 * classes to them. Arbitrarily large motif sizes are supported.
 *
 * \param graph The graph object to study.
 * \param no Pointer to an integer type, the result will be stored
 *        here.
 * \param size The size of the motifs to count.
 * \param cut_prob Vector giving the probabilities that a branch of
 *        the search tree will be cut at a given level.
 * \return Error code.
 * \sa \ref igraph_motifs_randesu(), \ref
 *     igraph_motifs_randesu_estimate().
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_motifs_randesu_no(const igraph_t *graph, igraph_integer_t *no,
                             igraph_integer_t size, const igraph_vector_t *cut_prob) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t neis;
    igraph_vector_int_t vids;     /* this is G */
    igraph_vector_int_t adjverts; /* this is V_E */
    igraph_stack_int_t stack;     /* this is S */
    igraph_integer_t *added;
    igraph_integer_t father;
    igraph_integer_t i;

    if (size < 3) {
        IGRAPH_ERRORF("Motif size must be at least 3, received %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, size);
    }

    if (igraph_vector_size(cut_prob) != size) {
        IGRAPH_ERRORF("Cut probability vector size (%" IGRAPH_PRId ") must agree with motif size (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL, igraph_vector_size(cut_prob), size);
    }
    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(added, "Insufficient memory to count motifs.");
    IGRAPH_FINALLY(igraph_free, added);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&adjverts, 0);
    IGRAPH_CHECK(igraph_stack_int_init(&stack, 0));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &stack);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    *no = 0;

    RNG_BEGIN();

    for (father = 0; father < no_of_nodes; father++) {
        igraph_integer_t level, s;

        IGRAPH_ALLOW_INTERRUPTION();

        if (VECTOR(*cut_prob)[0] == 1 ||
            RNG_UNIF01() < VECTOR(*cut_prob)[0]) {
            continue;
        }

        /* init G */
        igraph_vector_int_clear(&vids); level = 0;
        IGRAPH_CHECK(igraph_vector_int_push_back(&vids, father));
        added[father] += 1; level += 1;

        /* init V_E */
        igraph_vector_int_clear(&adjverts);
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
        s = igraph_vector_int_size(&neis);
        for (i = 0; i < s; i++) {
            igraph_integer_t nei = VECTOR(neis)[i];
            if (!added[nei] && nei > father) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei));
                IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, father));
            }
            added[nei] += 1;
        }

        /* init S */
        igraph_stack_int_clear(&stack);

        while (level > 1 || !igraph_vector_int_empty(&adjverts)) {
            igraph_real_t cp = VECTOR(*cut_prob)[level];

            if (level == size - 1) {
                s = igraph_vector_int_size(&adjverts) / 2;
                for (i = 0; i < s; i++) {
                    if (cp != 0 && RNG_UNIF01() < cp) {
                        continue;
                    }
                    (*no) += 1;
                }
            }

            if (level < size - 1 &&
                !igraph_vector_int_empty(&adjverts)) {
                /* We might step down */
                igraph_integer_t neifather = igraph_vector_int_pop_back(&adjverts);
                igraph_integer_t nei = igraph_vector_int_pop_back(&adjverts);

                if (cp == 0 || RNG_UNIF01() > cp) {
                    /* Yes, step down */
                    IGRAPH_CHECK(igraph_vector_int_push_back(&vids, nei));
                    added[nei] += 1; level += 1;

                    IGRAPH_CHECK(igraph_stack_int_push(&stack, neifather));
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, nei));
                    IGRAPH_CHECK(igraph_stack_int_push(&stack, level));

                    IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
                    s = igraph_vector_int_size(&neis);
                    for (i = 0; i < s; i++) {
                        igraph_integer_t nei2 = VECTOR(neis)[i];
                        if (!added[nei2] && nei2 > father) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei2));
                            IGRAPH_CHECK(igraph_vector_int_push_back(&adjverts, nei));
                        }
                        added[nei2] += 1;
                    }
                }
            } else {
                /* no, step back */
                igraph_integer_t nei, neifather;
                while (!igraph_stack_int_empty(&stack) &&
                       level == igraph_stack_int_top(&stack) - 1) {
                    igraph_stack_int_pop(&stack);
                    nei = igraph_stack_int_pop(&stack);
                    neifather = igraph_stack_int_pop(&stack);
                    igraph_vector_int_push_back(&adjverts, nei);
                    igraph_vector_int_push_back(&adjverts, neifather);
                }

                nei = igraph_vector_int_pop_back(&vids);
                added[nei] -= 1; level -= 1;
                IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
                s = igraph_vector_int_size(&neis);
                for (i = 0; i < s; i++) {
                    added[ VECTOR(neis)[i] ] -= 1;
                }
                while (!igraph_vector_int_empty(&adjverts) &&
                       igraph_vector_int_tail(&adjverts) == nei) {
                    igraph_vector_int_pop_back(&adjverts);
                    igraph_vector_int_pop_back(&adjverts);
                }
            }

        } /* while */

        /* clear the added vector */
        added[father] -= 1;
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
        s = igraph_vector_int_size(&neis);
        for (i = 0; i < s; i++) {
            added[ VECTOR(neis)[i] ] -= 1;
        }

    } /* for father */

    RNG_END();

    IGRAPH_FREE(added);
    igraph_vector_int_destroy(&vids);
    igraph_vector_int_destroy(&adjverts);
    igraph_stack_int_destroy(&stack);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(5);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_dyad_census
 * \brief Dyad census, as defined by Holland and Leinhardt.
 *
 * Dyad census means classifying each pair of vertices of a directed
 * graph into three categories: mutual (there is at least one edge from
 * \c a to \c b and also from \c b to \c a); asymmetric (there is at least
 * one edge either from \c a to \c b or from \c b to \c a, but not the other
 * way) and null (no edges between \c a and \c b in either direction).
 *
 * </para><para>
 * Holland, P.W. and Leinhardt, S.  (1970).  A Method for Detecting
 * Structure in Sociometric Data.  American Journal of Sociology,
 * 70, 492-513.
 *
 * \param graph The input graph. For an undirected graph, there are no
 *    asymmetric connections.
 * \param mut Pointer to a real, the number of mutual dyads is
 *    stored here.
 * \param asym Pointer to a real, the number of asymmetric dyads
 *    is stored here.
 * \param null Pointer to a real, the number of null dyads is
 *    stored here.
 * \return Error code.
 *
 * \sa \ref igraph_reciprocity(), \ref igraph_triad_census().
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */
igraph_error_t igraph_dyad_census(const igraph_t *graph, igraph_real_t *mut,
                       igraph_real_t *asym, igraph_real_t *null) {

    /* This function operates with a floating point type instead of an
     * integer type in order to avoid integer overflow, which is likely
     * for 'null' in large graphs on 32-bit systems. */

    igraph_real_t nonrec = 0, rec = 0;
    igraph_vector_int_t inneis, outneis;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&inneis, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&outneis, 0);

    for (i = 0; i < vc; i++) {
        igraph_integer_t ideg, odeg;
        igraph_integer_t ip, op;

        IGRAPH_CHECK(igraph_i_neighbors(graph, &inneis, i, IGRAPH_IN, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
        IGRAPH_CHECK(igraph_i_neighbors(graph, &outneis, i, IGRAPH_OUT, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));

        ideg = igraph_vector_int_size(&inneis);
        odeg = igraph_vector_int_size(&outneis);

        ip = op = 0;
        while (ip < ideg && op < odeg) {
            if (VECTOR(inneis)[ip] < VECTOR(outneis)[op]) {
                nonrec += 1;
                ip++;
            } else if (VECTOR(inneis)[ip] > VECTOR(outneis)[op]) {
                nonrec += 1;
                op++;
            } else {
                rec += 1;
                ip++;
                op++;
            }
        }
        nonrec += (ideg - ip) + (odeg - op);
    }

    igraph_vector_int_destroy(&inneis);
    igraph_vector_int_destroy(&outneis);
    IGRAPH_FINALLY_CLEAN(2);

    *mut = rec / 2;
    *asym = nonrec / 2;
    *null = 0.5 * vc * (vc - 1.0) - (*mut + *asym);
    if (*null == 0.0) *null = 0.0; /* avoid returning -0.0 */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_triad_census_24(const igraph_t *graph, igraph_real_t *res2,
                           igraph_real_t *res4) {

    igraph_integer_t vc = igraph_vcount(graph);
    igraph_vector_int_t seen;
    igraph_vector_int_t *neis, *neis2;
    igraph_integer_t i, j, k, s, neilen, neilen2, ign;
    igraph_adjlist_t adjlist;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&seen, vc);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    *res2 = *res4 = 0;

    for (i = 0; i < vc; i++) {
        IGRAPH_ALLOW_INTERRUPTION();

        neis = igraph_adjlist_get(&adjlist, i);
        neilen = igraph_vector_int_size(neis);
        /* mark neighbors of i & i itself */
        VECTOR(seen)[i] = i + 1;
        ign = 0;
        for (j = 0; j < neilen; j++) {
            igraph_integer_t nei = VECTOR(*neis)[j];
            if (VECTOR(seen)[nei] == i + 1 || VECTOR(seen)[nei] == -(i + 1)) {
                /* multiple edges or loop edge */
                VECTOR(seen)[nei] = -(i + 1);
                ign++;
            } else {
                VECTOR(seen)[nei] = i + 1;
            }
        }

        for (j = 0; j < neilen; j++) {
            igraph_integer_t nei = VECTOR(*neis)[j];
            if (nei <= i || (j > 0 && nei == VECTOR(*neis)[j - 1])) {
                continue;
            }
            neis2 = igraph_adjlist_get(&adjlist, nei);
            neilen2 = igraph_vector_int_size(neis2);
            s = 0;
            for (k = 0; k < neilen2; k++) {
                igraph_integer_t nei2 = VECTOR(*neis2)[k];
                if (k > 0 && nei2 == VECTOR(*neis2)[k - 1]) {
                    continue;
                }
                if (VECTOR(seen)[nei2] != i + 1 && VECTOR(seen)[nei2] != -(i + 1)) {
                    s++;
                }
            }
            if (VECTOR(seen)[nei] > 0) {
                *res2 += vc - s - neilen + ign - 1;
            } else {
                *res4 += vc - s - neilen + ign - 1;
            }
        }
    }

    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&seen);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_triad_census
 * \brief Triad census, as defined by Davis and Leinhardt.
 *
 * Calculating the triad census means classifying every triple of
 * vertices in a directed graph based on the type of pairwise
 * connections it contains, i.e. mutual, asymmetric or no connection.
 * A triple can be in one of 16 states, commonly described using
 * Davis and Leinhardt's "MAN labels". The \p res vector will
 * contain the counts of these in the following order:
 *
 * \clist
 * \cli &#xa0;0: 003
 *      A, B, C, the empty graph.
 * \cli &#xa0;1: 012
 *      A->B, C, a graph with a single directed edge.
 * \cli &#xa0;2: 102
 *      A&lt;->B, C, a graph with a mutual connection between two vertices.
 * \cli &#xa0;3: 021D
 *      A&lt;-B->C, the binary out-tree.
 * \cli &#xa0;4: 021U
 *      A->B&lt;-C, the binary in-tree.
 * \cli &#xa0;5: 021C
 *      A->B->C, the directed line.
 * \cli &#xa0;6: 111D
 *      A&lt;->B&lt;-C.
 * \cli &#xa0;7: 111U
 *      A&lt;->B->C.
 * \cli &#xa0;8: 030T
 *      A->B&lt;-C, A->C.
 * \cli &#xa0;9: 030C
 *      A&lt;-B&lt;-C, A->C.
 * \cli 10: 201
 *      A&lt;->B&lt;->C.
 * \cli 11: 120D
 *      A&lt;-B->C, A&lt;->C.
 * \cli 12: 120U
 *      A->B&lt;-C, A&lt;->C.
 * \cli 13: 120C
 *      A->B->C, A&lt;->C.
 * \cli 14: 210
 *      A->B&lt;->C, A&lt;->C.
 * \cli 15: 300
 *      A&lt;->B&lt;->C, A&lt;->C, the complete graph.
 * \endclist
 *
 * </para><para>
 * This function is intended for directed graphs. If the input is undirected,
 * a warning is shown, and undirected edges will be interpreted as mutual.
 *
 * </para><para>
 * This function calls \ref igraph_motifs_randesu() which is an
 * implementation of the FANMOD motif finder tool, see \ref
 * igraph_motifs_randesu() for details. Note that the order of the
 * triads is not the same for \ref igraph_triad_census() and \ref
 * igraph_motifs_randesu().
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Davis, J.A. and Leinhardt, S.  (1972).  The Structure of
 * Positive Interpersonal Relations in Small Groups.  In J. Berger
 * (Ed.), Sociological Theories in Progress, Volume 2, 218-251.
 * Boston: Houghton Mifflin.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the result is stored
 *   here in the same order as given in the list above. Note that this
 *   order is different than the one used by \ref igraph_motifs_randesu().
 * \return Error code.
 *
 * \sa \ref igraph_motifs_randesu(), \ref igraph_dyad_census().
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_triad_census(const igraph_t *graph, igraph_vector_t *res) {

    igraph_vector_t cut_prob;
    igraph_real_t m2, m4;
    igraph_vector_t tmp;
    igraph_integer_t vc = igraph_vcount(graph);
    igraph_real_t total;

    if (!igraph_is_directed(graph)) {
        IGRAPH_WARNING("Triad census called on an undirected graph. All connections will be treated as mutual.");
    }

    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&cut_prob, 3); /* all zeros */
    IGRAPH_CHECK(igraph_vector_resize(res, 16));
    igraph_vector_null(res);
    IGRAPH_CHECK(igraph_motifs_randesu(graph, &tmp, 3, &cut_prob));
    IGRAPH_CHECK(igraph_i_triad_census_24(graph, &m2, &m4));

    total = ((igraph_real_t)vc) * (vc - 1);
    total *= (vc - 2);
    total /= 6;

    /* Reorder */
    if (igraph_is_directed(graph)) {
        VECTOR(tmp)[0] = 0;
        VECTOR(tmp)[1] = m2;
        VECTOR(tmp)[3] = m4;
        VECTOR(tmp)[0] = total - igraph_vector_sum(&tmp);

        VECTOR(*res)[0] = VECTOR(tmp)[0];
        VECTOR(*res)[1] = VECTOR(tmp)[1];
        VECTOR(*res)[2] = VECTOR(tmp)[3];
        VECTOR(*res)[3] = VECTOR(tmp)[6];
        VECTOR(*res)[4] = VECTOR(tmp)[2];
        VECTOR(*res)[5] = VECTOR(tmp)[4];
        VECTOR(*res)[6] = VECTOR(tmp)[5];
        VECTOR(*res)[7] = VECTOR(tmp)[9];
        VECTOR(*res)[8] = VECTOR(tmp)[7];
        VECTOR(*res)[9] = VECTOR(tmp)[11];
        VECTOR(*res)[10] = VECTOR(tmp)[10];
        VECTOR(*res)[11] = VECTOR(tmp)[8];
        VECTOR(*res)[12] = VECTOR(tmp)[13];
        VECTOR(*res)[13] = VECTOR(tmp)[12];
        VECTOR(*res)[14] = VECTOR(tmp)[14];
        VECTOR(*res)[15] = VECTOR(tmp)[15];
    } else {
        VECTOR(tmp)[0] = 0;
        VECTOR(tmp)[1] = m2;
        VECTOR(tmp)[0] = total - igraph_vector_sum(&tmp);

        VECTOR(*res)[0] = VECTOR(tmp)[0];
        VECTOR(*res)[2] = VECTOR(tmp)[1];
        VECTOR(*res)[10] = VECTOR(tmp)[2];
        VECTOR(*res)[15] = VECTOR(tmp)[3];
    }

    igraph_vector_destroy(&cut_prob);
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
