/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_cliques.h"

#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_constants.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_stack.h"

#include "cliques/cliquer_internal.h"
#include "core/interruption.h"
#include "core/set.h"

#include <string.h>    /* memset */

static void igraph_i_cliques_free_res(igraph_vector_ptr_t *res) {
    long i, n;

    n = igraph_vector_ptr_size(res);
    for (i = 0; i < n; i++) {
        if (VECTOR(*res)[i] != 0) {
            igraph_vector_destroy(VECTOR(*res)[i]);
            igraph_free(VECTOR(*res)[i]);
        }
    }
    igraph_vector_ptr_clear(res);
}

static int igraph_i_find_k_cliques(
        const igraph_t *graph,
        long int size,
        const igraph_real_t *member_storage,
        igraph_real_t **new_member_storage,
        long int old_clique_count,
        long int *clique_count,
        igraph_vector_t *neis,
        igraph_bool_t independent_vertices) {

    long int j, k, l, m, n, new_member_storage_size;
    const igraph_real_t *c1, *c2;
    igraph_real_t v1, v2;
    igraph_bool_t ok;

    /* Allocate the storage */
    *new_member_storage = IGRAPH_REALLOC(*new_member_storage,
                                         (size_t) (size * old_clique_count),
                                         igraph_real_t);
    if (*new_member_storage == 0) {
        IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
    }
    new_member_storage_size = size * old_clique_count;
    IGRAPH_FINALLY(igraph_free, *new_member_storage);

    m = n = 0;

    /* Now consider all pairs of i-1-cliques and see if they can be merged */
    for (j = 0; j < old_clique_count; j++) {
        for (k = j + 1; k < old_clique_count; k++) {
            IGRAPH_ALLOW_INTERRUPTION();

            /* Since cliques are represented by their vertex indices in increasing
             * order, two cliques can be merged iff they have exactly the same
             * indices excluding one AND there is an edge between the two different
             * vertices */
            c1 = member_storage + j * (size - 1);
            c2 = member_storage + k * (size - 1);
            /* Find the longest prefixes of c1 and c2 that are equal */
            for (l = 0; l < size - 1 && c1[l] == c2[l]; l++) {
                (*new_member_storage)[m++] = c1[l];
            }
            /* Now, if l == size-1, the two vectors are totally equal.
            This is a bug */
            if (l == size - 1) {
                IGRAPH_WARNING("possible bug in igraph_cliques");
                m = n;
            } else {
                /* Assuming that j<k, c1[l] is always less than c2[l], since cliques
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
                ok = 1;
                for (; l < size - 1; l++) {
                    if (c1[l] == c2[l]) {
                        (*new_member_storage)[m++] = c1[l];
                        ok = 0;
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
                 * one place, so the whole clique is invalid. */
                if (l != size - 1) {
                    /* Step back in new_member_storage */
                    m = n;
                } else {
                    /* v1 and v2 are the two different vertices. Check for an edge
                     * if we are looking for cliques and check for the absence of an
                     * edge if we are looking for independent vertex sets */
                    IGRAPH_CHECK(igraph_neighbors(graph, neis, (igraph_integer_t) v1,
                                                  IGRAPH_ALL));
                    l = igraph_vector_search(neis, 0, v2, 0);
                    if ((l && !independent_vertices) || (!l && independent_vertices)) {
                        /* Found a new clique, step forward in new_member_storage */
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
                                                         igraph_real_t);
                    if (*new_member_storage == 0) {
                        IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
                    }
                    new_member_storage_size *= 2;
                    IGRAPH_FINALLY(igraph_free, *new_member_storage);
                }
            }
        }
    }

    /* Calculate how many cliques have we found */
    *clique_count = n / size;

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/* Internal function for calculating cliques or independent vertex sets.
 * They are practically the same except that the complementer of the graph
 * should be used in the latter case.
 */
static int igraph_i_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
                            igraph_integer_t min_size, igraph_integer_t max_size,
                            igraph_bool_t independent_vertices) {

    igraph_integer_t no_of_nodes;
    igraph_vector_t neis;
    igraph_real_t *member_storage = 0, *new_member_storage, *c1;
    long int i, j, k, clique_count, old_clique_count;

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("directionality of edges is ignored for directed graphs");
    }

    no_of_nodes = igraph_vcount(graph);

    if (min_size < 0) {
        min_size = 0;
    }
    if (max_size > no_of_nodes || max_size <= 0) {
        max_size = no_of_nodes;
    }

    igraph_vector_ptr_clear(res);

    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_FINALLY(igraph_i_cliques_free_res, res);

    /* Will be resized later, if needed. */
    member_storage = IGRAPH_CALLOC(1, igraph_real_t);
    if (member_storage == 0) {
        IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, member_storage);

    /* Find all 1-cliques: every vertex will be a clique */
    new_member_storage = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    if (new_member_storage == 0) {
        IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, new_member_storage);

    for (i = 0; i < no_of_nodes; i++) {
        new_member_storage[i] = i;
    }
    clique_count = no_of_nodes;
    old_clique_count = 0;

    /* Add size 1 cliques if requested */
    if (min_size <= 1) {
        IGRAPH_CHECK(igraph_vector_ptr_resize(res, no_of_nodes));
        igraph_vector_ptr_null(res);
        for (i = 0; i < no_of_nodes; i++) {
            igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
            if (p == 0) {
                IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
            }
            IGRAPH_FINALLY(igraph_free, p);
            IGRAPH_CHECK(igraph_vector_init(p, 1));
            VECTOR(*p)[0] = i;
            VECTOR(*res)[i] = p;
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    for (i = 2; i <= max_size && clique_count > 1; i++) {

        /* Here new_member_storage contains the cliques found in the previous
           iteration. Save this into member_storage, might be needed later  */

        c1 = member_storage;
        member_storage = new_member_storage;
        new_member_storage = c1;
        old_clique_count = clique_count;

        IGRAPH_ALLOW_INTERRUPTION();

        /* Calculate the cliques */

        IGRAPH_FINALLY_CLEAN(2);
        IGRAPH_CHECK(igraph_i_find_k_cliques(graph, i, member_storage,
                                             &new_member_storage,
                                             old_clique_count,
                                             &clique_count,
                                             &neis,
                                             independent_vertices));
        IGRAPH_FINALLY(igraph_free, member_storage);
        IGRAPH_FINALLY(igraph_free, new_member_storage);

        /* Add the cliques just found to the result if requested */
        if (i >= min_size && i <= max_size) {
            for (j = 0, k = 0; j < clique_count; j++, k += i) {
                igraph_vector_t *p = IGRAPH_CALLOC(1, igraph_vector_t);
                if (p == 0) {
                    IGRAPH_ERROR("cliques failed", IGRAPH_ENOMEM);
                }
                IGRAPH_FINALLY(igraph_free, p);
                IGRAPH_CHECK(igraph_vector_init_copy(p, &new_member_storage[k], i));
                IGRAPH_FINALLY(igraph_vector_destroy, p);
                IGRAPH_CHECK(igraph_vector_ptr_push_back(res, p));
                IGRAPH_FINALLY_CLEAN(2);
            }
        }

    } /* i <= max_size && clique_count != 0 */

    igraph_free(member_storage);
    igraph_free(new_member_storage);
    igraph_vector_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(4); /* 3 here, +1 is igraph_i_cliques_free_res */

    return 0;
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
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, i.e. \p res will contain pointers to \ref igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_largest_cliques() and \ref igraph_clique_number().
 *
 * Time complexity: Exponential
 *
 * \example examples/simple/igraph_cliques.c
 */
int igraph_cliques(const igraph_t *graph, igraph_vector_ptr_t *res,
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
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \return Error code.
 *
 * \sa \ref igraph_cliques() and \ref igraph_cliques_callback()
 *
 * Time complexity: Exponential
 *
 */
int igraph_clique_size_hist(const igraph_t *graph, igraph_vector_t *hist,
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
 * callback function as a pointer to an \ref igraph_vector_t.  Destroying and
 * freeing this vector is left up to the user.  Use \ref igraph_vector_destroy()
 * to destroy it first, then free it using \ref igraph_free().
 *
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * \param graph The input graph.
 * \param min_size Integer giving the minimum size of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param cliquehandler_fn Callback function to be called for each clique.
 * See also \ref igraph_clique_handler_t.
 * \param arg Extra argument to supply to \p cliquehandler_fn.
 * \return Error code.
 *
 * \sa \ref igraph_cliques()
 *
 * Time complexity: Exponential
 *
 */
int igraph_cliques_callback(const igraph_t *graph,
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
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * Only positive integer vertex weights are supported.
 *
 * \param graph The input graph.
 * \param vertex_weights A vector of vertex weights. The current implementation
 *   will truncate all weights to their integer parts.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, i.e. \p res will contain pointers to \ref igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param min_weight Integer giving the minimum weight of the cliques to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_weight Integer giving the maximum weight of the cliques to be
 *   returned. If negative or zero, no upper bound will be used.
 * \param maximal If true, only maximal cliques will be returned
 * \return Error code.
 *
 * \sa \ref igraph_cliques(), \ref igraph_maximal_cliques()
 *
 * Time complexity: Exponential
 *
 */
int igraph_weighted_cliques(const igraph_t *graph,
                            const igraph_vector_t *vertex_weights, igraph_vector_ptr_t *res,
                            igraph_real_t min_weight, igraph_real_t max_weight, igraph_bool_t maximal) {
    return igraph_i_weighted_cliques(graph, vertex_weights, res, min_weight, max_weight, maximal);
}


/**
 * \function igraph_largest_weighted_cliques
 * \brief Finds the largest weight clique(s) in a graph.
 *
 * </para><para>
 * Finds the clique(s) having the largest weight in the graph.
 *
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * Only positive integer vertex weights are supported.
 *
 * \param graph The input graph.
 * \param vertex_weights A vector of vertex weights. The current implementation
 *   will truncate all weights to their integer parts.
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, i.e. \p res will contain pointers to \ref igraph_vector_t
 *   objects which contain the indices of vertices involved in a clique.
 *   The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \return Error code.
 *
 * \sa \ref igraph_weighted_cliques(), \ref igraph_weighted_clique_number(), \ref igraph_largest_cliques()
 *
 * Time complexity: TODO
 */
int igraph_largest_weighted_cliques(const igraph_t *graph,
                                    const igraph_vector_t *vertex_weights, igraph_vector_ptr_t *res) {
    return igraph_i_largest_weighted_cliques(graph, vertex_weights, res);
}


/**
 * \function igraph_weighted_clique_number
 * \brief Finds the weight of the largest weight clique in the graph.
 *
 * </para><para>The current implementation of this function
 * uses version 1.21 of the Cliquer library by Sampo Niskanen and
 * Patric R. J. Östergård, http://users.aalto.fi/~pat/cliquer.html
 *
 * Only positive integer vertex weights are supported.
 *
 * \param graph The input graph.
 * \param vertex_weights A vector of vertex weights. The current implementation
 *   will truncate all weights to their integer parts.
 * \param res The largest weight will be returned to the \c igraph_real_t
 *   pointed to by this variable.
 * \return Error code.
 *
 * \sa \ref igraph_weighted_cliques(), \ref igraph_largest_weighted_cliques(), \ref igraph_clique_number()
 *
 * Time complexity: TODO
 *
 */
int igraph_weighted_clique_number(const igraph_t *graph,
                                  const igraph_vector_t *vertex_weights, igraph_real_t *res) {
    return igraph_i_weighted_clique_number(graph, vertex_weights, res);
}

typedef int(*igraph_i_maximal_clique_func_t)(const igraph_vector_t*, void*, igraph_bool_t*);
typedef struct {
    igraph_vector_ptr_t* result;
    igraph_integer_t min_size;
    igraph_integer_t max_size;
} igraph_i_maximal_clique_data_t;

static int igraph_i_maximal_cliques(const igraph_t *graph, igraph_i_maximal_clique_func_t func, void* data);

static int igraph_i_maximal_or_largest_cliques_or_indsets(
        const igraph_t *graph,
        igraph_vector_ptr_t *res,
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
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, i.e. \p res will contain pointers to \ref igraph_vector_t
 *   objects which contain the indices of vertices involved in an independent
 *   vertex set. The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \param min_size Integer giving the minimum size of the sets to be
 *   returned. If negative or zero, no lower bound will be used.
 * \param max_size Integer giving the maximum size of the sets to be
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
int igraph_independent_vertex_sets(const igraph_t *graph,
                                   igraph_vector_ptr_t *res,
                                   igraph_integer_t min_size,
                                   igraph_integer_t max_size) {
    return igraph_i_cliques(graph, res, min_size, max_size, 1);
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
 * \param res Pointer to a pointer vector, the result will be stored
 *     here. It will be resized as needed.
 * \return Error code.
 *
 * \sa \ref igraph_independent_vertex_sets(), \ref
 * igraph_maximal_independent_vertex_sets().
 *
 * Time complexity: TODO
 */

int igraph_largest_independent_vertex_sets(const igraph_t *graph,
        igraph_vector_ptr_t *res) {
    return igraph_i_maximal_or_largest_cliques_or_indsets(graph, res, 0, 1, 0);
}

typedef struct igraph_i_max_ind_vsets_data_t {
    igraph_integer_t matrix_size;
    igraph_adjlist_t adj_list;         /* Adjacency list of the graph */
    igraph_vector_t deg;                 /* Degrees of individual nodes */
    igraph_set_t* buckets;               /* Bucket array */
    /* The IS value for each node. Still to be explained :) */
    igraph_integer_t* IS;
    igraph_integer_t largest_set_size;   /* Size of the largest set encountered */
    igraph_bool_t keep_only_largest;     /* True if we keep only the largest sets */
} igraph_i_max_ind_vsets_data_t;

static int igraph_i_maximal_independent_vertex_sets_backtrack(
        const igraph_t *graph,
        igraph_vector_ptr_t *res,
        igraph_i_max_ind_vsets_data_t *clqdata,
        igraph_integer_t level) {
    long int v1, v2, v3, c, j, k;
    igraph_vector_int_t *neis1, *neis2;
    igraph_bool_t f;
    igraph_integer_t j1;
    long int it_state;

    IGRAPH_ALLOW_INTERRUPTION();

    if (level >= clqdata->matrix_size - 1) {
        igraph_integer_t size = 0;
        if (res) {
            igraph_vector_t *vec;
            vec = IGRAPH_CALLOC(1, igraph_vector_t);
            if (vec == 0) {
                IGRAPH_ERROR("igraph_i_maximal_independent_vertex_sets failed", IGRAPH_ENOMEM);
            }
            IGRAPH_VECTOR_INIT_FINALLY(vec, 0);
            for (v1 = 0; v1 < clqdata->matrix_size; v1++)
                if (clqdata->IS[v1] == 0) {
                    IGRAPH_CHECK(igraph_vector_push_back(vec, v1));
                }
            size = (igraph_integer_t) igraph_vector_size(vec);
            if (!clqdata->keep_only_largest) {
                IGRAPH_CHECK(igraph_vector_ptr_push_back(res, vec));
            } else {
                if (size > clqdata->largest_set_size) {
                    /* We are keeping only the largest sets, and we've found one that's
                     * larger than all previous sets, so we have to clear the list */
                    j = igraph_vector_ptr_size(res);
                    for (v1 = 0; v1 < j; v1++) {
                        igraph_vector_destroy(VECTOR(*res)[v1]);
                        free(VECTOR(*res)[v1]);
                    }
                    igraph_vector_ptr_clear(res);
                    IGRAPH_CHECK(igraph_vector_ptr_push_back(res, vec));
                } else if (size == clqdata->largest_set_size) {
                    IGRAPH_CHECK(igraph_vector_ptr_push_back(res, vec));
                } else {
                    igraph_vector_destroy(vec);
                    free(vec);
                }
            }
            IGRAPH_FINALLY_CLEAN(1);
        } else {
            for (v1 = 0, size = 0; v1 < clqdata->matrix_size; v1++)
                if (clqdata->IS[v1] == 0) {
                    size++;
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
               (v2 = (long int) VECTOR(*neis1)[j]) <= level) {
            if (clqdata->IS[v2] == 0) {
                c++;
            }
            j++;
        }

        if (c == 0) {
            /* If there are no such nodes... */
            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = (long int) VECTOR(*neis1)[j]) <= level) {
                clqdata->IS[v2]++;
                j++;
            }
            IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, clqdata, (igraph_integer_t) v1));
            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = (long int) VECTOR(*neis1)[j]) <= level) {
                clqdata->IS[v2]--;
                j++;
            }
        } else {
            /* If there are such nodes, store the count in the IS value of v1 */
            clqdata->IS[v1] = (igraph_integer_t) c;
            IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, clqdata, (igraph_integer_t) v1));
            clqdata->IS[v1] = 0;

            f = 1;
            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = (long int) VECTOR(*neis1)[j]) <= level) {
                if (clqdata->IS[v2] == 0) {
                    IGRAPH_CHECK(igraph_set_add(&clqdata->buckets[v1],
                                                (igraph_integer_t) j));
                    neis2 = igraph_adjlist_get(&clqdata->adj_list, v2);
                    k = 0;
                    while (k < VECTOR(clqdata->deg)[v2] &&
                           (v3 = (long int) VECTOR(*neis2)[k]) <= level) {
                        clqdata->IS[v3]--;
                        if (clqdata->IS[v3] == 0) {
                            f = 0;
                        }
                        k++;
                    }
                }
                clqdata->IS[v2]++;
                j++;
            }

            if (f) {
                IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, clqdata, (igraph_integer_t) v1));
            }

            j = 0;
            while (j < VECTOR(clqdata->deg)[v1] &&
                   (v2 = (long int) VECTOR(*neis1)[j]) <= level) {
                clqdata->IS[v2]--;
                j++;
            }

            it_state = 0;
            while (igraph_set_iterate(&clqdata->buckets[v1], &it_state, &j1)) {
                j = (long)j1;
                v2 = (long int) VECTOR(*neis1)[j];
                neis2 = igraph_adjlist_get(&clqdata->adj_list, v2);
                k = 0;
                while (k < VECTOR(clqdata->deg)[v2] &&
                       (v3 = (long int) VECTOR(*neis2)[k]) <= level) {
                    clqdata->IS[v3]++;
                    k++;
                }
            }
            igraph_set_clear(&clqdata->buckets[v1]);
        }
    }

    return 0;
}

static void igraph_i_free_set_array(igraph_set_t* array) {
    long int i = 0;
    while (igraph_set_inited(array + i)) {
        igraph_set_destroy(array + i);
        i++;
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
 * \param res Pointer to a pointer vector, the result will be stored
 *   here, i.e. \p res will contain pointers to \ref igraph_vector_t
 *   objects which contain the indices of vertices involved in an independent
 *   vertex set. The pointer vector will be resized if needed but note that the
 *   objects in the pointer vector will not be freed.
 * \return Error code.
 *
 * \sa \ref igraph_maximal_cliques(), \ref
 * igraph_independence_number()
 *
 * Time complexity: TODO.
 */
int igraph_maximal_independent_vertex_sets(const igraph_t *graph,
        igraph_vector_ptr_t *res) {
    igraph_i_max_ind_vsets_data_t clqdata;
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph), i;

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("directionality of edges is ignored for directed graphs");
    }

    clqdata.matrix_size = no_of_nodes;
    clqdata.keep_only_largest = 0;

    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &clqdata.adj_list, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &clqdata.adj_list);

    clqdata.IS = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (clqdata.IS == 0) {
        IGRAPH_ERROR("igraph_maximal_independent_vertex_sets failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, clqdata.IS);

    IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(clqdata.deg)[i] = igraph_vector_int_size(igraph_adjlist_get(&clqdata.adj_list, i));
    }

    clqdata.buckets = IGRAPH_CALLOC(no_of_nodes + 1, igraph_set_t);
    if (clqdata.buckets == 0) {
        IGRAPH_ERROR("igraph_maximal_independent_vertex_sets failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));
    }

    igraph_vector_ptr_clear(res);

    /* Do the show */
    clqdata.largest_set_size = 0;
    IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, &clqdata, 0));

    /* Cleanup */
    for (i = 0; i < no_of_nodes; i++) {
        igraph_set_destroy(&clqdata.buckets[i]);
    }
    igraph_adjlist_destroy(&clqdata.adj_list);
    igraph_vector_destroy(&clqdata.deg);
    igraph_free(clqdata.IS);
    igraph_free(clqdata.buckets);
    IGRAPH_FINALLY_CLEAN(4);
    return 0;
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
int igraph_independence_number(const igraph_t *graph, igraph_integer_t *no) {
    igraph_i_max_ind_vsets_data_t clqdata;
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph), i;

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("directionality of edges is ignored for directed graphs");
    }

    clqdata.matrix_size = no_of_nodes;
    clqdata.keep_only_largest = 0;

    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &clqdata.adj_list, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &clqdata.adj_list);

    clqdata.IS = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (clqdata.IS == 0) {
        IGRAPH_ERROR("igraph_independence_number failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, clqdata.IS);

    IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(clqdata.deg)[i] = igraph_vector_int_size(igraph_adjlist_get(&clqdata.adj_list, i));
    }

    clqdata.buckets = IGRAPH_CALLOC(no_of_nodes + 1, igraph_set_t);
    if (clqdata.buckets == 0) {
        IGRAPH_ERROR("igraph_independence_number failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));
    }

    /* Do the show */
    clqdata.largest_set_size = 0;
    IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, 0, &clqdata, 0));
    *no = clqdata.largest_set_size;

    /* Cleanup */
    for (i = 0; i < no_of_nodes; i++) {
        igraph_set_destroy(&clqdata.buckets[i]);
    }
    igraph_adjlist_destroy(&clqdata.adj_list);
    igraph_vector_destroy(&clqdata.deg);
    igraph_free(clqdata.IS);
    igraph_free(clqdata.buckets);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

/*************************************************************************/
/* MAXIMAL CLIQUES, LARGEST CLIQUES                                      */
/*************************************************************************/

static int igraph_i_maximal_cliques_store_max_size(const igraph_vector_t* clique, void* data,
                                                   igraph_bool_t* cont) {
    igraph_integer_t* result = (igraph_integer_t*)data;
    IGRAPH_UNUSED(cont);
    if (*result < igraph_vector_size(clique)) {
        *result = (igraph_integer_t) igraph_vector_size(clique);
    }
    return IGRAPH_SUCCESS;
}

/*
static int igraph_i_maximal_cliques_store(const igraph_vector_t* clique, void* data, igraph_bool_t* cont) {
    igraph_vector_ptr_t* result = (igraph_vector_ptr_t*)data;
    igraph_vector_t* vec;

    IGRAPH_UNUSED(cont);
    vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (vec == 0) {
        IGRAPH_ERROR("cannot allocate memory for storing next clique", IGRAPH_ENOMEM);
    }

    IGRAPH_CHECK(igraph_vector_copy(vec, clique));
    IGRAPH_CHECK(igraph_vector_ptr_push_back(result, vec));

    return IGRAPH_SUCCESS;
}

static int igraph_i_maximal_cliques_store_size_check(const igraph_vector_t* clique, void* data_, igraph_bool_t* cont) {
    igraph_i_maximal_clique_data_t* data = (igraph_i_maximal_clique_data_t*)data_;
    igraph_vector_t* vec;
    igraph_integer_t size = (igraph_integer_t) igraph_vector_size(clique);

    IGRAPH_UNUSED(cont);
    if (size < data->min_size || size > data->max_size) {
        return IGRAPH_SUCCESS;
    }

    vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (vec == 0) {
        IGRAPH_ERROR("cannot allocate memory for storing next clique", IGRAPH_ENOMEM);
    }

    IGRAPH_CHECK(igraph_vector_copy(vec, clique));
    IGRAPH_CHECK(igraph_vector_ptr_push_back(data->result, vec));

    return IGRAPH_SUCCESS;
}
*/

static int igraph_i_largest_cliques_store(const igraph_vector_t* clique, void* data, igraph_bool_t* cont) {
    igraph_vector_ptr_t* result = (igraph_vector_ptr_t*)data;
    igraph_vector_t* vec;
    long int i, n;

    IGRAPH_UNUSED(cont);
    /* Is the current clique at least as large as the others that we have found? */
    if (!igraph_vector_ptr_empty(result)) {
        n = igraph_vector_size(clique);
        if (n < igraph_vector_size(VECTOR(*result)[0])) {
            return IGRAPH_SUCCESS;
        }

        if (n > igraph_vector_size(VECTOR(*result)[0])) {
            for (i = 0; i < igraph_vector_ptr_size(result); i++) {
                igraph_vector_destroy(VECTOR(*result)[i]);
            }
            igraph_vector_ptr_free_all(result);
            igraph_vector_ptr_resize(result, 0);
        }
    }

    vec = IGRAPH_CALLOC(1, igraph_vector_t);
    if (vec == 0) {
        IGRAPH_ERROR("cannot allocate memory for storing next clique", IGRAPH_ENOMEM);
    }

    IGRAPH_CHECK(igraph_vector_copy(vec, clique));
    IGRAPH_CHECK(igraph_vector_ptr_push_back(result, vec));

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
 * for maximal cliques using \ref igraph_maximal_cliques() and drops
 * those that are not the largest.
 *
 * </para><para>The implementation of this function changed between
 * igraph 0.5 and 0.6, so the order of the cliques and the order of
 * vertices within the cliques will almost surely be different between
 * these two versions.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized pointer vector, the result
 *        will be stored here. It will be resized as needed. Note that
 *        vertices of a clique may be returned in arbitrary order.
 * \return Error code.
 *
 * \sa \ref igraph_cliques(), \ref igraph_maximal_cliques()
 *
 * Time complexity: O(3^(|V|/3)) worst case.
 */

int igraph_largest_cliques(const igraph_t *graph, igraph_vector_ptr_t *res) {
    igraph_vector_ptr_clear(res);
    IGRAPH_FINALLY(igraph_i_cliques_free_res, res);
    IGRAPH_CHECK(igraph_i_maximal_cliques(graph, &igraph_i_largest_cliques_store, (void*)res));
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_clique_number
 * \brief Finds the clique number of the graph.
 *
 * </para><para>
 * The clique number of a graph is the size of the largest clique.
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
int igraph_clique_number(const igraph_t *graph, igraph_integer_t *no) {
    *no = 0;
    return igraph_i_maximal_cliques(graph, &igraph_i_maximal_cliques_store_max_size, (void*)no);
}

typedef struct {
    igraph_vector_int_t cand;
    igraph_vector_int_t fini;
    igraph_vector_int_t cand_filtered;
} igraph_i_maximal_cliques_stack_frame;

static void igraph_i_maximal_cliques_stack_frame_destroy(igraph_i_maximal_cliques_stack_frame *frame) {
    igraph_vector_int_destroy(&frame->cand);
    igraph_vector_int_destroy(&frame->fini);
    igraph_vector_int_destroy(&frame->cand_filtered);
}

static void igraph_i_maximal_cliques_stack_destroy(igraph_stack_ptr_t *stack) {
    igraph_i_maximal_cliques_stack_frame *frame;

    while (!igraph_stack_ptr_empty(stack)) {
        frame = (igraph_i_maximal_cliques_stack_frame*)igraph_stack_ptr_pop(stack);
        igraph_i_maximal_cliques_stack_frame_destroy(frame);
        free(frame);
    }

    igraph_stack_ptr_destroy(stack);
}

static int igraph_i_maximal_cliques(const igraph_t *graph, igraph_i_maximal_clique_func_t func, void* data) {
    int directed = igraph_is_directed(graph);
    long int i, j, k, l;
    igraph_integer_t no_of_nodes, nodes_to_check, nodes_done;
    igraph_integer_t best_cand = 0, best_cand_degree = 0, best_fini_cand_degree;
    igraph_adjlist_t adj_list;
    igraph_stack_ptr_t stack;
    igraph_i_maximal_cliques_stack_frame frame, *new_frame_ptr;
    igraph_vector_t clique;
    igraph_vector_int_t new_cand, new_fini, cn, best_cand_nbrs,
                        best_fini_cand_nbrs;
    igraph_bool_t cont = 1;
    igraph_bool_t found;

    if (directed) {
        IGRAPH_WARNING("directionality of edges is ignored for directed graphs");
    }

    no_of_nodes = igraph_vcount(graph);
    if (no_of_nodes == 0) {
        return IGRAPH_SUCCESS;
    }

    /* Construct an adjacency list representation */
    IGRAPH_CHECK(igraph_adjlist_init(
        graph, &adj_list, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adj_list);
    igraph_adjlist_sort(&adj_list);

    /* Initialize stack */
    IGRAPH_CHECK(igraph_stack_ptr_init(&stack, 0));
    IGRAPH_FINALLY(igraph_i_maximal_cliques_stack_destroy, &stack);

    /* Create the initial (empty) clique */
    IGRAPH_VECTOR_INIT_FINALLY(&clique, 0);

    /* Initialize new_cand, new_fini, cn, best_cand_nbrs and best_fini_cand_nbrs (will be used later) */
    igraph_vector_int_init(&new_cand, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &new_cand);
    igraph_vector_int_init(&new_fini, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &new_fini);
    igraph_vector_int_init(&cn, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &cn);
    igraph_vector_int_init(&best_cand_nbrs, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &best_cand_nbrs);
    igraph_vector_int_init(&best_fini_cand_nbrs, 0);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &best_fini_cand_nbrs);

    /* Find the vertex with the highest degree */
    best_cand = 0; best_cand_degree = (igraph_integer_t) igraph_vector_int_size(igraph_adjlist_get(&adj_list, 0));
    for (i = 1; i < no_of_nodes; i++) {
        j = igraph_vector_int_size(igraph_adjlist_get(&adj_list, i));
        if (j > best_cand_degree) {
            best_cand = (igraph_integer_t) i;
            best_cand_degree = (igraph_integer_t) j;
        }
    }

    /* Create the initial stack frame */
    IGRAPH_CHECK(igraph_vector_int_init_seq(&frame.cand, 0, no_of_nodes - 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &frame.cand);
    IGRAPH_CHECK(igraph_vector_int_init(&frame.fini, 0));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &frame.fini);
    IGRAPH_CHECK(igraph_vector_int_init(&frame.cand_filtered, 0));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &frame.cand_filtered);
    IGRAPH_CHECK(igraph_vector_int_difference_sorted(&frame.cand,
                 igraph_adjlist_get(&adj_list, best_cand), &frame.cand_filtered));
    IGRAPH_FINALLY_CLEAN(3);
    IGRAPH_FINALLY(igraph_i_maximal_cliques_stack_frame_destroy, &frame);

    /* TODO: frame.cand and frame.fini should be a set instead of a vector */

    /* Main loop starts here */
    nodes_to_check = (igraph_integer_t) igraph_vector_int_size(&frame.cand_filtered); nodes_done = 0;
    while (!igraph_vector_int_empty(&frame.cand_filtered) || !igraph_stack_ptr_empty(&stack)) {
        if (igraph_vector_int_empty(&frame.cand_filtered)) {
            /* No candidates left to check in this stack frame, pop out the previous stack frame */
            igraph_i_maximal_cliques_stack_frame *newframe = igraph_stack_ptr_pop(&stack);
            igraph_i_maximal_cliques_stack_frame_destroy(&frame);
            frame = *newframe;
            free(newframe);

            if (igraph_stack_ptr_size(&stack) == 1) {
                /* We will be using the next candidate node in the next iteration, so we can increase
                 * nodes_done by 1 */
                nodes_done++;
            }

            /* For efficiency reasons, we only check for interruption and show progress here */
            IGRAPH_PROGRESS("Maximal cliques: ", 100.0 * nodes_done / nodes_to_check, NULL);
            IGRAPH_ALLOW_INTERRUPTION();

            igraph_vector_pop_back(&clique);
            continue;
        }

        /* Try the next node in the clique */
        i = (long int) igraph_vector_int_pop_back(&frame.cand_filtered);
        IGRAPH_CHECK(igraph_vector_push_back(&clique, i));

        /* Remove the node from the candidate list */
        found = igraph_vector_int_binsearch(&frame.cand, i, &j); IGRAPH_ASSERT(found);
        igraph_vector_int_remove(&frame.cand, j);

        /* Add the node to the finished list */
        found = igraph_vector_int_binsearch(&frame.fini, i, &j); IGRAPH_ASSERT(!found);
        IGRAPH_CHECK(igraph_vector_int_insert(&frame.fini, j, i));

        /* Create new_cand and new_fini */
        IGRAPH_CHECK(igraph_vector_int_intersect_sorted(&frame.cand, igraph_adjlist_get(&adj_list, i), &new_cand));
        IGRAPH_CHECK(igraph_vector_int_intersect_sorted(&frame.fini, igraph_adjlist_get(&adj_list, i), &new_fini));

        /* Do we have anything more to search? */
        if (igraph_vector_int_empty(&new_cand)) {
            if (igraph_vector_int_empty(&new_fini)) {
                /* We have a maximal clique here */
                IGRAPH_CHECK(func(&clique, data, &cont));
                if (!cont) {
                    /* The callback function requested to stop the search */
                    break;
                }
            }
            igraph_vector_pop_back(&clique);
            continue;
        }
        if (igraph_vector_int_empty(&new_fini) &&
            igraph_vector_int_size(&new_cand) == 1) {
            /* Shortcut: only one node left */
            IGRAPH_CHECK(igraph_vector_push_back(&clique, VECTOR(new_cand)[0]));
            IGRAPH_CHECK(func(&clique, data, &cont));
            if (!cont) {
                /* The callback function requested to stop the search */
                break;
            }
            igraph_vector_pop_back(&clique);
            igraph_vector_pop_back(&clique);
            continue;
        }

        /* Find the next best candidate node in new_fini */
        l = igraph_vector_int_size(&new_cand);
        best_cand_degree = -1;
        j = igraph_vector_int_size(&new_fini);
        for (i = 0; i < j; i++) {
            k = (long int)VECTOR(new_fini)[i];
            IGRAPH_CHECK(igraph_vector_int_intersect_sorted(&new_cand, igraph_adjlist_get(&adj_list, k), &cn));
            if (igraph_vector_int_size(&cn) > best_cand_degree) {
                best_cand_degree = (igraph_integer_t) igraph_vector_int_size(&cn);
                IGRAPH_CHECK(igraph_vector_int_update(&best_fini_cand_nbrs, &cn));
                if (best_cand_degree == l) {
                    /* Cool, we surely have the best candidate node here as best_cand_degree can't get any better */
                    break;
                }
            }
        }
        /* Shortcut here: we don't have to examine new_cand */
        if (best_cand_degree == l) {
            igraph_vector_pop_back(&clique);
            continue;
        }
        /* Still finding best candidate node */
        best_fini_cand_degree = best_cand_degree;
        best_cand_degree = -1;
        j = igraph_vector_int_size(&new_cand);
        l = l - 1;
        for (i = 0; i < j; i++) {
            k = (long int)VECTOR(new_cand)[i];
            IGRAPH_CHECK(igraph_vector_int_intersect_sorted(&new_cand, igraph_adjlist_get(&adj_list, k), &cn));
            if (igraph_vector_int_size(&cn) > best_cand_degree) {
                best_cand_degree = (igraph_integer_t) igraph_vector_int_size(&cn);
                IGRAPH_CHECK(igraph_vector_int_update(&best_cand_nbrs, &cn));
                if (best_cand_degree == l) {
                    /* Cool, we surely have the best candidate node here as best_cand_degree can't get any better */
                    break;
                }
            }
        }

        /* Create a new stack frame in case we back out later */
        new_frame_ptr = IGRAPH_CALLOC(1, igraph_i_maximal_cliques_stack_frame);
        if (new_frame_ptr == 0) {
            IGRAPH_ERROR("cannot allocate new stack frame", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, new_frame_ptr);
        *new_frame_ptr = frame;
        memset(&frame, 0, sizeof(frame));
        IGRAPH_CHECK(igraph_stack_ptr_push(&stack, new_frame_ptr));
        IGRAPH_FINALLY_CLEAN(1);  /* ownership of new_frame_ptr taken by the stack */
        /* Ownership of the current frame and its vectors (frame.cand, frame.done, frame.cand_filtered)
         * is taken by the stack from now on. Vectors in frame must be re-initialized with new_cand,
         * new_fini and stuff. The old frame.cand and frame.fini won't be leaked because they are
         * managed by the stack now. */
        frame.cand = new_cand;
        frame.fini = new_fini;
        IGRAPH_CHECK(igraph_vector_int_init(&new_cand, 0));
        IGRAPH_CHECK(igraph_vector_int_init(&new_fini, 0));
        IGRAPH_CHECK(igraph_vector_int_init(&frame.cand_filtered, 0));

        /* Adjust frame.cand_filtered */
        if (best_cand_degree < best_fini_cand_degree) {
            IGRAPH_CHECK(igraph_vector_int_difference_sorted(&frame.cand, &best_fini_cand_nbrs, &frame.cand_filtered));
        } else {
            IGRAPH_CHECK(igraph_vector_int_difference_sorted(&frame.cand, &best_cand_nbrs, &frame.cand_filtered));
        }
    }

    IGRAPH_PROGRESS("Maximal cliques: ", 100.0, NULL);

    igraph_adjlist_destroy(&adj_list);
    igraph_vector_destroy(&clique);
    igraph_vector_int_destroy(&new_cand);
    igraph_vector_int_destroy(&new_fini);
    igraph_vector_int_destroy(&cn);
    igraph_vector_int_destroy(&best_cand_nbrs);
    igraph_vector_int_destroy(&best_fini_cand_nbrs);
    igraph_i_maximal_cliques_stack_frame_destroy(&frame);
    igraph_i_maximal_cliques_stack_destroy(&stack);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

static int igraph_i_maximal_or_largest_cliques_or_indsets(const igraph_t *graph,
        igraph_vector_ptr_t *res,
        igraph_integer_t *clique_number,
        igraph_bool_t keep_only_largest,
        igraph_bool_t complementer) {
    igraph_i_max_ind_vsets_data_t clqdata;
    igraph_integer_t no_of_nodes = (igraph_integer_t) igraph_vcount(graph), i;

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("directionality of edges is ignored for directed graphs");
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
    if (clqdata.IS == 0) {
        IGRAPH_ERROR("igraph_i_maximal_or_largest_cliques_or_indsets failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, clqdata.IS);

    IGRAPH_VECTOR_INIT_FINALLY(&clqdata.deg, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(clqdata.deg)[i] = igraph_vector_int_size(igraph_adjlist_get(&clqdata.adj_list, i));
    }

    clqdata.buckets = IGRAPH_CALLOC(no_of_nodes + 1, igraph_set_t);
    if (clqdata.buckets == 0) {
        IGRAPH_ERROR("igraph_maximal_or_largest_cliques_or_indsets failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_free_set_array, clqdata.buckets);

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_CHECK(igraph_set_init(&clqdata.buckets[i], 0));
    }

    if (res) {
        igraph_vector_ptr_clear(res);
    }

    /* Do the show */
    clqdata.largest_set_size = 0;
    IGRAPH_CHECK(igraph_i_maximal_independent_vertex_sets_backtrack(graph, res, &clqdata, 0));

    /* Cleanup */
    for (i = 0; i < no_of_nodes; i++) {
        igraph_set_destroy(&clqdata.buckets[i]);
    }
    igraph_adjlist_destroy(&clqdata.adj_list);
    igraph_vector_destroy(&clqdata.deg);
    igraph_free(clqdata.IS);
    igraph_free(clqdata.buckets);
    IGRAPH_FINALLY_CLEAN(4);

    if (clique_number) {
        *clique_number = clqdata.largest_set_size;
    }
    return 0;
}
