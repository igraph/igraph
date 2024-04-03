/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2023 The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_cliques.h"

#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"

#include "cliques/cliquer_internal.h"
#include "core/interruption.h"
#include "core/set.h"

static igraph_error_t igraph_i_find_k_indsets(
        const igraph_t *graph,
        igraph_integer_t size,
        const igraph_integer_t *member_storage,
        igraph_integer_t **new_member_storage,
        igraph_integer_t old_count,
        igraph_integer_t *new_count,
        igraph_vector_int_t *neis) {

    igraph_integer_t l, m, n, new_member_storage_size;
    const igraph_integer_t *c1, *c2;
    igraph_integer_t v1, v2;
    igraph_bool_t ok;

    /* Allocate the storage */
    *new_member_storage = IGRAPH_REALLOC(*new_member_storage,
                                         (size_t) (size * old_count),
                                         igraph_integer_t);
    IGRAPH_CHECK_OOM(*new_member_storage, "Insufficient memory for independent vertex sets.");

    new_member_storage_size = size * old_count;
    IGRAPH_FINALLY(igraph_free, *new_member_storage);

    m = n = 0;

    /* Now consider all pairs of i-1-indsets and see if they can be merged */
    for (igraph_integer_t j = 0; j < old_count; j++) {
        for (igraph_integer_t k = j + 1; k < old_count; k++) {
            IGRAPH_ALLOW_INTERRUPTION();

            /* Since indsets are represented by their vertex indices in increasing
             * order, two indsets can be merged iff they have exactly the same
             * indices excluding one AND there is no edge between the two different
             * vertices */
            c1 = member_storage + j * (size - 1);
            c2 = member_storage + k * (size - 1);
            /* Find the longest prefixes of c1 and c2 that are equal */
            for (l = 0; l < size - 1 && c1[l] == c2[l]; l++) {
                (*new_member_storage)[m++] = c1[l];
            }
            /* Now, if l == size-1, the two vectors are totally equal. This is a bug */
            IGRAPH_ASSERT(l != size-1);
            /* Assuming that j<k, c1[l] is always less than c2[l], since indsets
                * are ordered alphabetically. Now add c1[l] and store c2[l] in a
                * dummy variable */
            (*new_member_storage)[m++] = c1[l];
            v1 = c1[l];
            v2 = c2[l];
            l++;
            /* Copy the remaining part of the two vectors. Every member pair
                * found in the remaining parts satisfies the following:
                * 1. If they are equal, they should be added.
                * 2. If they are not equal, the smaller must be equal to the
                *    one stored in the dummy variable. If not, the two vectors
                *    differ in more than one place. The larger will be stored in
                *    the dummy variable again.
                */
            ok = true;
            for (; l < size - 1; l++) {
                if (c1[l] == c2[l]) {
                    (*new_member_storage)[m++] = c1[l];
                    ok = false;
                } else if (ok) {
                    if (c1[l] < c2[l]) {
                        if (c1[l] == v1) {
                            (*new_member_storage)[m++] = c1[l];
                            v2 = c2[l];
                        } else {
                            break;
                        }
                    } else {
                        if (ok && c2[l] == v1) {
                            (*new_member_storage)[m++] = c2[l];
                            v2 = c1[l];
                        } else {
                            break;
                        }
                    }
                } else {
                    break;
                }
            }
            /* Now, if l != size-1, the two vectors had a difference in more than
                * one place, so the whole independent vertex set is invalid. */
            if (l != size - 1) {
                /* Step back in new_member_storage */
                m = n;
            } else {
                /* v1 and v2 are the two different vertices. Check for the
                    * absence of an edge since we are looking for independent
                    * vertex sets */
                IGRAPH_CHECK(igraph_neighbors(graph, neis, v1, IGRAPH_ALL));
                if (!igraph_vector_int_search(neis, 0, v2, 0)) {
                    /* Found a new independent vertex set, step forward in new_member_storage */
                    if (m == n || v2 > (*new_member_storage)[m - 1]) {
                        (*new_member_storage)[m++] = v2;
                        n = m;
                    } else {
                        m = n;
                    }
                } else {
                    m = n;
                }
            }
            /* See if new_member_storage is full. If so, reallocate */
            if (m == new_member_storage_size) {
                IGRAPH_FINALLY_CLEAN(1);
                *new_member_storage = IGRAPH_REALLOC(*new_member_storage,
                                                        (size_t) new_member_storage_size * 2,
                                                        igraph_integer_t);
                IGRAPH_CHECK_OOM(*new_member_storage, "igraph_independent_vertex_sets failed");
                new_member_storage_size *= 2;
                IGRAPH_FINALLY(igraph_free, *new_member_storage);
            }
        }
    }

    /* Calculate how many independent vertex sets we have found */
    *new_count = n / size;

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_cliques
 * \brief Finds all or some cliques in a graph.
 *
 * </para><para>
 * Cliques are fully connected subgraphs of a graph.
 *
 * </para><para>
 * If you are only interested in the size of the largest clique in the graph,
 * use \ref igraph_clique_number() instead.
 *
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \param min_size Integer specifying the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer specifying the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_largest_cliques() and \ref igraph_clique_number().
 *
 * Time complexity: Exponential
 *
 * \example examples/simple/igraph_cliques.c
 */
igraph_error_t igraph_cliques(const igraph_t *graph, igraph_vector_int_list_t *res,
                   igraph_integer_t min_size, igraph_integer_t max_size) {
    return igraph_i_cliquer_cliques(graph, res, min_size, max_size);
}


/**
 * \function igraph_clique_size_hist
 * \brief Counts cliques of each size in the graph.
 *
 * </para><para>
 * Cliques are fully connected subgraphs of a graph.
 *
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param hist Pointer to an initialized vector. The result will be stored
 * here. The first element will store the number of size-1 cliques, the second
 * element the number of size-2 cliques, etc.  For cliques smaller than \p min_size,
 * zero counts will be returned.
 * \param min_size Integer specifying the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer specifying the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_cliques() and \ref igraph_cliques_callback()
 *
 * Time complexity: Exponential
 *
 */
igraph_error_t igraph_clique_size_hist(const igraph_t *graph, igraph_vector_t *hist,
                            igraph_integer_t min_size, igraph_integer_t max_size) {
    return igraph_i_cliquer_histogram(graph, hist, min_size, max_size);
}


/**
 * \function igraph_cliques_callback
 * \brief Calls a function for each clique in the graph.
 *
 * </para><para>
 * Cliques are fully connected subgraphs of a graph. This function
 * enumerates all cliques within the given size range and calls
 * \p cliquehandler_fn for each of them. The cliques are passed to the
 * callback function as a pointer to an \ref igraph_vector_int_t.  Destroying and
 * freeing this vector is left up to the user.  Use \ref igraph_vector_int_destroy()
 * to destroy it first, then free it using \ref igraph_free().
 *
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param min_size Integer specifying the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer specifying the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param cliquehandler_fn Callback function to be called for each clique.
 *   See also \ref igraph_clique_handler_t.
 * \param arg Extra argument to supply to \p cliquehandler_fn.
 * \return Error code.
 *
 * \sa \ref igraph_cliques()
 *
 * Time complexity: Exponential
 *
 */
igraph_error_t igraph_cliques_callback(const igraph_t *graph,
                            igraph_integer_t min_size, igraph_integer_t max_size,
                            igraph_clique_handler_t *cliquehandler_fn, void *arg) {
    return igraph_i_cliquer_callback(graph, min_size, max_size, cliquehandler_fn, arg);
}


/**
 * \function igraph_weighted_cliques
 * \brief Finds all cliques in a given weight range in a vertex weighted graph.
 *
 * </para><para>
 * Cliques are fully connected subgraphs of a graph.
 * The weight of a clique is the sum of the weights
 * of individual vertices within the clique.
 *
 * </para><para>
 * Only positive integer vertex weights are supported.
 *
 * </para><para>
 * The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param vertex_weights A vector of vertex weights. The current implementation
 *   will truncate all weights to their integer parts. You may pass \c NULL
 *   here to make each vertex have a weight of 1.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \param min_weight Integer specifying the minimum weight of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_weight Integer specifying the maximum weight of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param maximal If true, only maximal cliques will be returned
 * \return Error code.
 *
 * \sa \ref igraph_cliques(), \ref igraph_maximal_cliques()
 *
 * Time complexity: Exponential
 *
 */
igraph_error_t igraph_weighted_cliques(const igraph_t *graph,
                            const igraph_vector_t *vertex_weights, igraph_vector_int_list_t *res,
                            igraph_real_t min_weight, igraph_real_t max_weight, igraph_bool_t maximal) {
    if (vertex_weights) {
        return igraph_i_weighted_cliques(graph, vertex_weights, res, min_weight, max_weight, maximal);
    } else if (maximal) {
        return igraph_maximal_cliques(graph, res, min_weight, max_weight);
    } else {
        return igraph_cliques(graph, res, min_weight, max_weight);
    }
}


/**
 * \function igraph_largest_weighted_cliques
 * \brief Finds the largest weight clique(s) in a graph.
 *
 * The weight of a clique is the sum of the weights of its vertices.
 * This function finds the clique(s) having the largest weight in the graph.
 *
 * </para><para>
 * Only positive integer vertex weights are supported.
 *
 * </para><para>
 * The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param vertex_weights A vector of vertex weights. The current implementation
 *   will truncate all weights to their integer parts. You may pass \c NULL
 *   here to make each vertex have a weight of 1.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \return Error code.
 *
 * \sa \ref igraph_weighted_cliques(), \ref igraph_weighted_clique_number(), \ref igraph_largest_cliques()
 *
 * Time complexity: TODO
 */
igraph_error_t igraph_largest_weighted_cliques(const igraph_t *graph,
                                    const igraph_vector_t *vertex_weights, igraph_vector_int_list_t *res) {
    if (vertex_weights) {
        return igraph_i_largest_weighted_cliques(graph, vertex_weights, res);
    } else {
        return igraph_largest_cliques(graph, res);
    }
}


/**
 * \function igraph_weighted_clique_number
 * \brief Finds the weight of the largest weight clique in the graph.
 *
 * The weight of a clique is the sum of the weights of its vertices.
 * This function finds the weight of the largest weight clique.
 *
 * </para><para>
 * Only positive integer vertex weights are supported.
 *
 * </para><para>
 * The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param vertex_weights A vector of vertex weights. The current implementation
 *   will truncate all weights to their integer parts. You may pass \c NULL
 *   here to make each vertex have a weight of 1.
 * \param res The largest weight will be returned to the \c igraph_real_t
 *   pointed to by this variable.
 * \return Error code.
 *
 * \sa \ref igraph_weighted_cliques(), \ref igraph_largest_weighted_cliques(), \ref igraph_clique_number()
 *
 * Time complexity: TODO
 *
 */
igraph_error_t igraph_weighted_clique_number(const igraph_t *graph,
                                  const igraph_vector_t *vertex_weights, igraph_real_t *res) {
    if (vertex_weights) {
        return igraph_i_weighted_clique_number(graph, vertex_weights, res);
    } else {
        igraph_integer_t res_int;
        IGRAPH_CHECK(igraph_clique_number(graph, &res_int));
        if (res) {
            *res = res_int;
        }
        return IGRAPH_SUCCESS;
    }
}

static igraph_error_t igraph_i_maximal_or_largest_cliques_or_indsets(
        const igraph_t *graph,
        igraph_vector_int_list_t *res,
        igraph_integer_t *clique_number,
        igraph_bool_t keep_only_largest,
        igraph_bool_t complementer);

/**
 * \function igraph_independent_vertex_sets
 * \brief Finds all independent vertex sets in a graph.
 *
 * </para><para>
 * A vertex set is considered independent if there are no edges between
 * them.
 *
 * </para><para>
 * If you are interested in the size of the largest independent vertex set,
 * use \ref igraph_independence_number() instead.
 *
 * </para><para>
 * The current implementation was ported to igraph from the Very Nauty Graph
 * Library by Keith Briggs and uses the algorithm from the paper
 * S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm
 * for generating all the maximal independent sets. SIAM J Computing,
 * 6:505--517, 1977.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \param min_size Integer specifying the minimum size of the sets to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer specifying the maximum size of the sets to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_largest_independent_vertex_sets(),
 * \ref igraph_independence_number().
 *
 * Time complexity: TODO
 *
 * \example examples/simple/igraph_independent_sets.c
 */
igraph_error_t igraph_independent_vertex_sets(const igraph_t *graph,
                                   igraph_vector_int_list_t *res,
                                   igraph_integer_t min_size,
                                   igraph_integer_t max_size) {
    igraph_integer_t no_of_nodes;
    igraph_vector_int_t neis, *indset;
    igraph_integer_t *member_storage, *new_member_storage, *c1;
    igraph_vector_int_t new_member_storage_view;
    igraph_integer_t indset_count, old_indset_count;

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored during independent vertex set calculations.");
    }

    no_of_nodes = igraph_vcount(graph);

    if (min_size < 0) {
        min_size = 0;
    }
    if (max_size > no_of_nodes || max_size <= 0) {
        max_size = no_of_nodes;
    }

    igraph_vector_int_list_clear(res);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    /* Will be resized later, if needed. */
    member_storage = IGRAPH_CALLOC(1, igraph_integer_t);
    IGRAPH_CHECK_OOM(member_storage, "Insufficient memory for independent vertex set calculation.");
    IGRAPH_FINALLY(igraph_free, member_storage);

    /* Find all 1-cliques: every vertex will be a clique */
    new_member_storage = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(new_member_storage, "Insufficient memory for independent vertex set calculation.");
    IGRAPH_FINALLY(igraph_free, new_member_storage);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        new_member_storage[i] = i;
    }
    indset_count = no_of_nodes;
    old_indset_count = 0;

    /* Add size 1 indsets if requested */
    if (min_size <= 1) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(res, no_of_nodes));
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            indset = igraph_vector_int_list_get_ptr(res, i);
            IGRAPH_CHECK(igraph_vector_int_push_back(indset, i));
        }
    }

    for (igraph_integer_t i = 2; i <= max_size && indset_count > 1; i++) {

        /* Here new_member_storage contains the independent vertex sets found in
           the previous iteration. Save this into member_storage, might be needed later  */

        c1 = member_storage;
        member_storage = new_member_storage;
        new_member_storage = c1;
        old_indset_count = indset_count;

        IGRAPH_ALLOW_INTERRUPTION();

        /* Calculate the independent vertex sets */

        IGRAPH_FINALLY_CLEAN(2);
        IGRAPH_CHECK(igraph_i_find_k_indsets(graph, i, member_storage,
                                             &new_member_storage,
                                             old_indset_count,
                                             &indset_count,
                                             &neis));
        IGRAPH_FINALLY(igraph_free, member_storage);
        IGRAPH_FINALLY(igraph_free, new_member_storage);

        /* Add the cliques just found to the result if requested */
        if (i >= min_size && i <= max_size) {
            for (igraph_integer_t j = 0, k = 0; j < indset_count; j++, k += i) {
                igraph_vector_int_view(&new_member_storage_view, new_member_storage + k, i);
                IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, &new_member_storage_view));
            }
        }

    } /* i <= max_size && clique_count != 0 */

    IGRAPH_FREE(new_member_storage);
    IGRAPH_FREE(member_storage);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_largest_independent_vertex_sets
 * \brief Finds the largest independent vertex set(s) in a graph.
 *
 * </para><para>
 * An independent vertex set is largest if there is no other
 * independent vertex set with more vertices in the graph.
 *
 * </para><para>
 * The current implementation was ported to igraph from the Very Nauty Graph
 * Library by Keith Briggs and uses the algorithm from the paper
 * S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm
 * for generating all the maximal independent sets. SIAM J Computing,
 * 6:505--517, 1977.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \return Error code.
 *
 * \sa \ref igraph_independent_vertex_sets(), \ref
 * igraph_maximal_independent_vertex_sets().
 *
 * Time complexity: TODO
 */

igraph_error_t igraph_largest_independent_vertex_sets(const igraph_t *graph,
        igraph_vector_int_list_t *res) {
    return igraph_i_maximal_or_largest_cliques_or_indsets(graph, res, 0, true, false);
}

typedef struct igraph_i_max_ind_vsets_data_t {
    igraph_integer_t matrix_size;
    igraph_adjlist_t adj_list;         /* Adjacency list of the graph */
    igraph_vector_int_t deg;                 /* Degrees of individual nodes */
    igraph_set_t* buckets;               /* Bucket array */
    /* The IS value for each node. Still to be explained :) */
    igraph_integer_t* IS;
    igraph_integer_t largest_set_size;   /* Size of the largest set encountered */
    igraph_bool_t keep_only_largest;     /* True if we keep only the largest sets */
} igraph_i_max_ind_vsets_data_t;

static igraph_error_t igraph_i_maximal_independent_vertex_sets_backtrack(
        const igraph_t *graph,
        igraph_vector_int_list_t *res,
        igraph_i_max_ind_vsets_data_t *clqdata,
        igraph_integer_t level) {

    igraph_integer_t v1, v2, v3, c, j, k;
    igraph_vector_int_t *neis1, *neis2;
    igraph_bool_t f;
    igraph_integer_t it_state;

    IGRAPH_ALLOW_INTERRUPTION();

    if (level >= clqdata->matrix_size - 1) {
        igraph_integer_t size = 0;

        if (res) {
            igraph_vector_int_t vec, *newvec;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&vec, 0);

            for (v1 = 0; v1 < clqdata->matrix_size; v1++) {
                if (clqdata->IS[v1] == 0) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&vec, v1));
                }
            }

            size = igraph_vector_int_size(&vec);

            /* Trick for efficient insertion of a new vector into a vector list:
             * Instead of copying the vector contents, we add an empty vector to
             * the list, then swap it with the vector to-be-added in O(1) time. */
            if (!clqdata->keep_only_largest) {
                IGRAPH_CHECK(igraph_vector_int_list_push_back_new(res, &newvec));
                igraph_vector_int_swap(newvec, &vec);
            } else {
                if (size > clqdata->largest_set_size) {
                    /* We are keeping only the largest sets, and we've found one that's
                     * larger than all previous sets, so we have to clear the list */
                    igraph_vector_int_list_clear(res);
                    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(res, &newvec));
                    igraph_vector_int_swap(newvec, &vec);
                } else if (size == clqdata->largest_set_size) {
                    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(res, &newvec));
                    igraph_vector_int_swap(newvec, &vec);
                }
            }

            igraph_vector_int_destroy(&vec);
            IGRAPH_FINALLY_CLEAN(1);
        } else {
            for (v1 = 0, size = 0; v1 < clqdata->matrix_size; v1++) {
                if (clqdata->IS[v1] == 0) {
                    size++;
                }
            }
        }
        if (size > clqdata->largest_set_size) {
            clqdata->largest_set_size = size;
        }
    } else {
        v1 = level + 1;
        /* Count the number of vertices with an index less than v1 that have
         * an IS value of zero */
        neis1 = igraph_adjlist_get(&clqdata->adj_list, v1);
        c = 0;
        j = 0;
        while (j < VECTOR(clqdata->deg)[v1] &&
               (v2 = VECTOR(*neis1)[j]) <= level) {
            if (clqdata->IS[v2] == 0) {
                c++;
            }
            j++;
        }

        if (c == 0) {
            /* If there are no such nodes... */
            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = VECTOR(*neis1)[j]) <= level) {
                clqdata->IS[v2]++;
                j++;
            }
            IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, clqdata, v1));
            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = VECTOR(*neis1)[j]) <= level) {
                clqdata->IS[v2]--;
                j++;
            }
        } else {
            /* If there are such nodes, store the count in the IS value of v1 */
            clqdata->IS[v1] = c;
            IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, clqdata, v1));
            clqdata->IS[v1] = 0;

            f = true;
            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = VECTOR(*neis1)[j]) <= level) {
                if (clqdata->IS[v2] == 0) {
                    IGRAPH_CHECK(igraph_set_add(&clqdata->buckets[v1], j));
                    neis2 = igraph_adjlist_get(&clqdata->adj_list, v2);
                    k = 0;
                    while (k < VECTOR(clqdata->deg)[v2] &&
                           (v3 = VECTOR(*neis2)[k]) <= level) {
                        clqdata->IS[v3]--;
                        if (clqdata->IS[v3] == 0) {
                            f = false;
                        }
                        k++;
                    }
                }
                clqdata->IS[v2]++;
                j++;
            }

            if (f) {
                IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, clqdata, v1));
            }

            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = VECTOR(*neis1)[j]) <= level) {
                clqdata->IS[v2]--;
                j++;
            }

            it_state = 0;
            while (igraph_set_iterate(&clqdata->buckets[v1], &it_state, &j)) {
                v2 = VECTOR(*neis1)[j];
                neis2 = igraph_adjlist_get(&clqdata->adj_list, v2);
                k = 0;
                while (k < VECTOR(clqdata->deg)[v2] &&
                       (v3 = VECTOR(*neis2)[k]) <= level) {
                    clqdata->IS[v3]++;
                    k++;
                }
            }
            igraph_set_clear(&clqdata->buckets[v1]);
        }
    }

    return IGRAPH_SUCCESS;
}

/* TODO (ugly hack):
 *
 * This version does not know the length of the array, and is safe to use
 * ONLY on arrays which have not been completely filled out and were
 * originally initialized to zero. It relies on igraph_set_inited()
 * returning false when igraph_set_t is all-zero-bytes.
 * This function is meant for use with IGRAPH_FINALLY.
 *
 * Should probably be replaced with a proper igraph_vector_ptr_t.
 */
static void free_set_array_incomplete(igraph_set_t *array) {
    igraph_integer_t i = 0;
    while (igraph_set_inited(array + i)) {
        igraph_set_destroy(array + i);
        i++;
    }
    IGRAPH_FREE(array);
}

static void free_set_array(igraph_set_t *array, igraph_integer_t n) {
    for (igraph_integer_t i=0; i < n; i++) {
        igraph_set_destroy(&array[i]);
    }
    IGRAPH_FREE(array);
}

/**
 * \function igraph_maximal_independent_vertex_sets
 * \brief Finds all maximal independent vertex sets of a graph.
 *
 * </para><para>
 * A maximal independent vertex set is an independent vertex set which
 * can't be extended any more by adding a new vertex to it.
 *
 * </para><para>
 * The algorithm used here is based on the following paper:
 * S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm for
 * generating all the maximal independent sets. SIAM J Computing,
 * 6:505--517, 1977.
 *
 * </para><para>
 * The implementation was originally written by Kevin O'Neill and modified
 * by K M Briggs in the Very Nauty Graph Library. I simply re-wrote it to
 * use igraph's data structures.
 *
 * </para><para>
 * If you are interested in the size of the largest independent vertex set,
 * use \ref igraph_independence_number() instead.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques(), \ref
 * igraph_independence_number()
 *
 * Time complexity: TODO.
 */
igraph_error_t igraph_maximal_independent_vertex_sets(const igraph_t *graph,
        igraph_vector_int_list_t *res) {
    igraph_i_max_ind_vsets_data_t clqdata;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored during independent vertex set calculations.");
    }

    clqdata.matrix_size = no_of_nodes;
    clqdata.keep_only_largest = false;

    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &clqdata.adj_list, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &clqdata.adj_list);

    clqdata.IS = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(clqdata.IS, "Insufficient memory for maximal independent vertex sets.");
    IGRAPH_FINALLY(igraph_free, clqdata.IS);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&clqdata.deg, no_of_nodes);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        VECTOR(clqdata.deg)[i] = igraph_vector_int_size(igraph_adjlist_get(&clqdata.adj_list, i));
    }

    clqdata.buckets = IGRAPH_CALLOC(no_of_nodes + 1, igraph_set_t);
    IGRAPH_CHECK_OOM(clqdata.buckets, "Insufficient memory for maximal independent vertex sets.");
    IGRAPH_FINALLY(free_set_array_incomplete, clqdata.buckets);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));
    }

    igraph_vector_int_list_clear(res);

    /* Do the show */
    clqdata.largest_set_size = 0;
    IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, &clqdata, 0));

    /* Cleanup */
    free_set_array(clqdata.buckets, no_of_nodes);
    igraph_vector_int_destroy(&clqdata.deg);
    IGRAPH_FREE(clqdata.IS);
    igraph_adjlist_destroy(&clqdata.adj_list);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_independence_number
 * \brief Finds the independence number of the graph.
 *
 * </para><para>
 * The independence number of a graph is the cardinality of the largest
 * independent vertex set.
 *
 * </para><para>
 * The current implementation was ported to igraph from the Very Nauty Graph
 * Library by Keith Briggs and uses the algorithm from the paper
 * S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm
 * for generating all the maximal independent sets. SIAM J Computing,
 * 6:505--517, 1977.
 *
 * \param graph The input graph.
 * \param no The independence number will be returned to the \c
 *   igraph_integer_t pointed by this variable.
 * \return Error code.
 *
 * \sa \ref igraph_independent_vertex_sets().
 *
 * Time complexity: TODO.
 */
igraph_error_t igraph_independence_number(const igraph_t *graph, igraph_integer_t *no) {
    igraph_i_max_ind_vsets_data_t clqdata;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored during independence number calculations.");
    }

    clqdata.matrix_size = no_of_nodes;
    clqdata.keep_only_largest = false;

    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &clqdata.adj_list, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &clqdata.adj_list);

    clqdata.IS = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(clqdata.IS, "Insufficient memory for independence number calculation.");
    IGRAPH_FINALLY(igraph_free, clqdata.IS);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&clqdata.deg, no_of_nodes);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        VECTOR(clqdata.deg)[i] = igraph_vector_int_size(igraph_adjlist_get(&clqdata.adj_list, i));
    }

    clqdata.buckets = IGRAPH_CALLOC(no_of_nodes + 1, igraph_set_t);
    IGRAPH_CHECK_OOM(clqdata.buckets, "Insufficient memory for independence number calculation.");
    IGRAPH_FINALLY(free_set_array_incomplete, clqdata.buckets);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));
    }

    /* Do the show */
    clqdata.largest_set_size = 0;
    IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, 0, &clqdata, 0));
    *no = clqdata.largest_set_size;

    /* Cleanup */
    free_set_array(clqdata.buckets, no_of_nodes);
    igraph_vector_int_destroy(&clqdata.deg);
    IGRAPH_FREE(clqdata.IS);
    igraph_adjlist_destroy(&clqdata.adj_list);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/*************************************************************************/
/* MAXIMAL CLIQUES, LARGEST CLIQUES                                      */
/*************************************************************************/

static igraph_error_t igraph_i_maximal_cliques_store_max_size(const igraph_vector_int_t* clique, void* data) {
    igraph_integer_t* result = (igraph_integer_t*)data;
    if (*result < igraph_vector_int_size(clique)) {
        *result = igraph_vector_int_size(clique);
    }
    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_largest_cliques_store(const igraph_vector_int_t* clique, void* data) {
    igraph_vector_int_list_t* result = (igraph_vector_int_list_t*)data;
    igraph_integer_t n;

    /* Is the current clique at least as large as the others that we have found? */
    if (!igraph_vector_int_list_empty(result)) {
        igraph_vector_int_t* first;

        n = igraph_vector_int_size(clique);
        first = igraph_vector_int_list_get_ptr(result, 0);
        if (n < igraph_vector_int_size(first)) {
            return IGRAPH_SUCCESS;
        }

        if (n > igraph_vector_int_size(first)) {
            igraph_vector_int_list_clear(result);
        }
    }

    IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(result, clique));

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_largest_cliques
 * \brief Finds the largest clique(s) in a graph.
 *
 * </para><para>
 * A clique is largest (quite intuitively) if there is no other clique
 * in the graph which contains more vertices.
 *
 * </para><para>
 * Note that this is not necessarily the same as a maximal clique,
 * i.e. the largest cliques are always maximal but a maximal clique is
 * not always largest.
 *
 * </para><para>The current implementation of this function searches
 * for maximal cliques using \ref igraph_maximal_cliques_callback() and drops
 * those that are not the largest.
 *
 * </para><para>The implementation of this function changed between
 * igraph 0.5 and 0.6, so the order of the cliques and the order of
 * vertices within the cliques will almost surely be different between
 * these two versions.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized list of integer vectors. The cliques
 *   will be stored here as vectors of vertex IDs.
 * \return Error code.
 *
 * \sa \ref igraph_cliques(), \ref igraph_maximal_cliques()
 *
 * Time complexity: O(3^(|V|/3)) worst case.
 */

igraph_error_t igraph_largest_cliques(const igraph_t *graph, igraph_vector_int_list_t *res) {
    igraph_vector_int_list_clear(res);
    IGRAPH_CHECK(igraph_maximal_cliques_callback(graph, &igraph_i_largest_cliques_store, (void*)res, 0, 0));
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_clique_number
 * \brief Finds the clique number of the graph.
 *
 * </para><para>
 * The clique number of a graph is the size of the largest clique.
 *
 * </para><para>The current implementation of this function searches
 * for maximal cliques using \ref igraph_maximal_cliques_callback() and keeps
 * track of the size of the largest clique that was found.
 *
 * \param graph The input graph.
 * \param no The clique number will be returned to the \c igraph_integer_t
 *   pointed by this variable.
 * \return Error code.
 *
 * \sa \ref igraph_cliques(), \ref igraph_largest_cliques().
 *
 * Time complexity: O(3^(|V|/3)) worst case.
 */
igraph_error_t igraph_clique_number(const igraph_t *graph, igraph_integer_t *no) {
    *no = 0;
    return igraph_maximal_cliques_callback(graph, &igraph_i_maximal_cliques_store_max_size, (void*)no, 0, 0);
}

static igraph_error_t igraph_i_maximal_or_largest_cliques_or_indsets(const igraph_t *graph,
        igraph_vector_int_list_t *res,
        igraph_integer_t *clique_number,
        igraph_bool_t keep_only_largest,
        igraph_bool_t complementer) {

    igraph_i_max_ind_vsets_data_t clqdata;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored for largest independent vertex set or clique calculations.");
    }

    clqdata.matrix_size = no_of_nodes;
    clqdata.keep_only_largest = keep_only_largest;

    if (complementer) {
        IGRAPH_CHECK(igraph_adjlist_init_complementer(graph, &clqdata.adj_list, IGRAPH_ALL, 0));
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(
            graph, &clqdata.adj_list, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE
        ));
    }
    IGRAPH_FINALLY(igraph_adjlist_destroy, &clqdata.adj_list);

    clqdata.IS = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(clqdata.IS, "Insufficient memory for largest independent sets or cliques.");
    IGRAPH_FINALLY(igraph_free, clqdata.IS);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&clqdata.deg, no_of_nodes);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        VECTOR(clqdata.deg)[i] = igraph_vector_int_size(igraph_adjlist_get(&clqdata.adj_list, i));
    }

    clqdata.buckets = IGRAPH_CALLOC(no_of_nodes + 1, igraph_set_t);
    IGRAPH_CHECK_OOM(clqdata.buckets, "Insufficient memory for largest independent sets or cliques.");
    IGRAPH_FINALLY(free_set_array_incomplete, clqdata.buckets);

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));
    }

    if (res) {
        igraph_vector_int_list_clear(res);
    }

    /* Do the show */
    clqdata.largest_set_size = 0;
    IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, &clqdata, 0));

    /* Cleanup */
    free_set_array(clqdata.buckets, no_of_nodes);
    igraph_vector_int_destroy(&clqdata.deg);
    igraph_free(clqdata.IS);
    igraph_adjlist_destroy(&clqdata.adj_list);
    IGRAPH_FINALLY_CLEAN(4);

    if (clique_number) {
        *clique_number = clqdata.largest_set_size;
    }

    return IGRAPH_SUCCESS;
}
